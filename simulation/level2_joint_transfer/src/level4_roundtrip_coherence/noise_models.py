"""可选技术噪声：quasistatic 与 Ornstein-Uhlenbeck 分数强度噪声。

所有 RMS、相关时间与相关系数均来自配置并标记 assumed（敏感性专用），
不声称由论文确定。噪声使用独立 seed，与初态 seed 分离。
"""
from __future__ import annotations

import numpy as np


class IntensityNoise:
    """每 shot 的 SLM/AOD 分数深度噪声（quasistatic 或 OU）。

    quasistatic：ε ~ N(0, σ²) 每 shot 一次抽样，SLM/AOD 可相关。
    OU：ε(t+dt) = ε·exp(-dt/τ) + sqrt(1-exp(-2dt/τ))·σ·n，边界处从
    quasistatic 抽样热启动。
    """

    def __init__(self, model, n_shots, rms, tau_s=None, cross_correlation=0.0,
                 seed=0):
        self.model = model
        self.n_shots = int(n_shots)
        self.rms = float(rms)
        self.tau_s = None if tau_s is None else float(tau_s)
        self.rho = float(np.clip(cross_correlation, -1.0, 1.0))
        self.rng = np.random.default_rng(int(seed))
        self._draw_initial()

    def _correlated_normals(self, n):
        """相关系数 rho 的两束独立高斯。"""
        z1 = self.rng.normal(size=n)
        z2 = self.rng.normal(size=n)
        if abs(self.rho) < 1e-12:
            return z1, z2
        return z1, self.rho * z1 + np.sqrt(1.0 - self.rho ** 2) * z2

    def _draw_initial(self):
        z1, z2 = self._correlated_normals(self.n_shots)
        self.eps_slm = self.rms * z1
        self.eps_aod = self.rms * z2

    def advance(self, dt):
        """OU 更新一步；quasistatic 恒定。"""
        if self.model != "ornstein_uhlenbeck" or self.tau_s is None:
            return
        decay = np.exp(-dt / self.tau_s)
        z1, z2 = self._correlated_normals(self.n_shots)
        scale = np.sqrt(max(0.0, 1.0 - decay ** 2)) * self.rms
        self.eps_slm = self.eps_slm * decay + scale * z1
        self.eps_aod = self.eps_aod * decay + scale * z2

    def depth_multipliers(self):
        """返回 (1+ε_slm, 1+ε_aod) 数组。"""
        return 1.0 + self.eps_slm, 1.0 + self.eps_aod

    def describe(self):
        """噪声参数描述（含 assumed 标注）。"""
        return {
            "model": self.model, "rms": self.rms,
            "tau_s": self.tau_s, "cross_correlation": self.rho,
            "provenance": "assumed_sensitivity_only",
        }


def noise_statistics(noise: IntensityNoise, n_samples=200000):
    """抽样验证噪声的均值/方差/相关时间（测试用）。"""
    probe = IntensityNoise(noise.model, n_samples, noise.rms, noise.tau_s,
                           noise.rho, seed=noise.rng.integers(1 << 31))
    eps0 = probe.eps_slm.copy()
    if noise.model == "ornstein_uhlenbeck" and noise.tau_s:
        dt = noise.tau_s / 50.0
        probe.advance(dt)
        corr = float(np.corrcoef(eps0, probe.eps_slm)[0, 1])
        tau_measured = -dt / np.log(corr) if 0 < corr < 1 else None
    else:
        tau_measured = None
    return {
        "mean_slm": float(np.mean(probe.eps_slm)),
        "std_slm": float(np.std(probe.eps_slm)),
        "mean_aod": float(np.mean(probe.eps_aod)),
        "std_aod": float(np.std(probe.eps_aod)),
        "cross_corr": float(np.corrcoef(probe.eps_slm, probe.eps_aod)[0, 1]),
        "tau_measured_s": tau_measured,
    }
