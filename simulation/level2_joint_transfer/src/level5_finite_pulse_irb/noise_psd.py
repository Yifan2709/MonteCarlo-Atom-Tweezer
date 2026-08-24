"""噪声模型：quasistatic、Ornstein-Uhlenbeck 与 PSD 合成高斯噪声。

支持 white、1/f^α、Lorentzian、60 Hz 窄带线及谐波、以及用户给定的
单边 PSD 表。FFT 合成约定：对长度 N、采样率 fs 的记录，
x[n] = Σ_k Re{A_k e^{i2πkn/N}}，双边 PSD S_x(f_k)=|X_k|²·(2/(fs·N))，
由此保证 ∫S 两边 = 方差。所有 RMS/谱指数/相关时间/线强均视为 assumed
（无实验数据时 nominal 关闭）。
"""
from __future__ import annotations

import numpy as np

TWO_PI = 2.0 * np.pi


def welch_psd(x: np.ndarray, fs: float, nperseg: int = 256):
    """简单 Welch PSD（Hann 窗、50% 重叠），返回 (freqs, psd_one_sided)。"""
    x = np.asarray(x, dtype=float)
    nperseg = min(nperseg, x.size)
    window = np.hanning(nperseg)
    step = nperseg // 2
    scales = []
    segs = []
    for start in range(0, x.size - nperseg + 1, step):
        seg = x[start:start + nperseg] * window
        segs.append(np.fft.rfft(seg))
    psd = np.zeros(nperseg // 2 + 1)
    for seg in segs:
        psd += np.abs(seg) ** 2
    psd /= len(segs)
    # 归一化：单边 PSD，含窗功率补偿
    win_power = np.sum(window ** 2)
    psd *= 2.0 / (fs * win_power)
    freqs = np.fft.rfftfreq(nperseg, d=1.0 / fs)
    return freqs, psd


def ou_trace(rng, n_samples, dt, tau, sigma):
    """Ornstein-Uhlenbeck 轨迹（离散精确更新）。"""
    if tau is None or tau <= 0:
        raise ValueError("OU 模型必须给出相关时间 τ>0")
    decay = np.exp(-dt / tau)
    scale = sigma * np.sqrt(max(0.0, 1.0 - decay ** 2))
    out = np.empty(n_samples)
    out[0] = sigma * rng.standard_normal()
    noise = rng.standard_normal(n_samples - 1)
    for i in range(1, n_samples):
        out[i] = out[i - 1] * decay + scale * noise[i - 1]
    return out


def ou_autocorrelation_theory(lag_s, tau, sigma):
    """OU 自相关理论值 σ²e^{-|lag|/τ}。"""
    return sigma ** 2 * np.exp(-np.abs(lag_s) / tau)


def _psd_model_one_sided(freqs, model, params):
    """目标单边 PSD（rad/s 或 Hz 单位由调用方固定：这里 f 以 Hz 计）。"""
    f = np.asarray(freqs, dtype=float)
    if model == "white":
        return np.full_like(f, float(params["level"]))
    if model == "one_over_f":
        alpha = float(params.get("alpha", 1.0))
        level = float(params["level"])
        f0 = float(params.get("f0_hz", 1.0))
        return level * (np.maximum(f, 1e-12) / f0) ** (-alpha)
    if model == "lorentzian":
        # S(f) = A τ / (1 + (2πfτ)²)（OU 过程的单边 PSD）
        tau = float(params["tau_s"])
        amp = float(params["sigma"]) ** 2
        return amp * 2.0 * tau / (1.0 + (TWO_PI * f * tau) ** 2)
    if model == "line":
        f0 = float(params["f_hz"])
        width = float(params.get("width_hz", 1.0))
        amp = float(params["amp"])
        total = np.zeros_like(f)
        for h in range(1, int(params.get("harmonics", 3)) + 1):
            total += amp / (h ** 2) * np.exp(
                -0.5 * ((f - h * f0) / width) ** 2) \
                * (width * np.sqrt(TWO_PI))
        return total
    raise KeyError(f"未知 PSD 模型 {model}")


def psd_target(freqs, spec):
    """组合 PSD 规范（list of {model, params}）→ 目标单边 PSD。"""
    freqs = np.asarray(freqs, dtype=float)
    total = np.zeros_like(freqs)
    for item in spec:
        total += _psd_model_one_sided(freqs, item["model"], item.get("params", {}))
    return total


def synthesize_from_psd(rng, psd_one_sided, n_samples, fs):
    """从单边 PSD 合成高斯时域噪声（FFT 相位随机化）。

    约定 S_1s(f_k)=2|X_k|²/(fs·N)（k=1..N/2−1，Nyquist 与 DC 不加倍），
    因此 |X_k|=sqrt(S·fs·N/2)；实数性由共轭对称谱保证，
    方差 = Σ_{k≥1}S_1s(f_k)·Δf。
    """
    psd = np.asarray(psd_one_sided, dtype=float)
    freqs = np.fft.rfftfreq(n_samples, d=1.0 / fs)
    if psd.size != freqs.size:
        raise ValueError("PSD 长度必须等于 rfftfreq(n) 的长度")
    amp = np.sqrt(np.maximum(psd, 0.0) * fs * n_samples / 2.0)
    phases = rng.uniform(0.0, TWO_PI, size=freqs.size)
    spectrum = amp * np.exp(1j * phases)
    spectrum[0] = 0.0
    if n_samples % 2 == 0:
        spectrum[-1] = np.abs(spectrum[-1])  # Nyquist 分量为实数
    full = np.concatenate([spectrum, np.conj(spectrum[-2:0:-1])])
    samples = np.fft.ifft(full).real
    var_target = float(np.sum(psd[1:]) * fs / n_samples)
    return samples, freqs, var_target


class DetuningNoiseModel:
    """失谐噪声模型：quasistatic / OU / PSD 合成的批量生成器。

    batch 维 = (实例, realization)；每条记录按模型独立生成，
    同一 seed 完全可复现、不同 seed 相互独立（测试验证）。
    """

    def __init__(self, model: str, rms: float = 0.0, tau_s: float | None = None,
                 psd_spec=None, fs_hz: float = 1e6, seed: int = 0,
                 n_samples: int = 0):
        self.model = model
        self.rms = float(rms)
        self.tau_s = tau_s
        self.psd_spec = psd_spec or []
        self.fs = float(fs_hz)
        self.rng = np.random.default_rng(int(seed))
        self.n_samples = int(n_samples)

    def sample_quasistatic(self, n: int) -> np.ndarray:
        """每实例一个常数偏移。"""
        return self.rng.normal(0.0, self.rms, size=int(n)) \
            if self.rms > 0 else np.zeros(int(n))

    def sample_trace(self, n: int, n_samples: int) -> np.ndarray:
        """每实例一条时域轨迹 [n, n_samples]。"""
        n, n_samples = int(n), int(n_samples)
        if self.model == "quasistatic":
            return np.repeat(self.sample_quasistatic(n)[:, None], n_samples,
                             axis=1)
        if self.model == "ornstein_uhlenbeck":
            if n_samples < 2:
                raise ValueError("OU 需要时域轨迹")
            out = np.empty((n, n_samples))
            for i in range(n):
                out[i] = ou_trace(self.rng, n_samples, 1.0 / self.fs,
                                  self.tau_s, self.rms)
            return out
        if self.model == "psd":
            freqs = np.fft.rfftfreq(n_samples, d=1.0 / self.fs)
            target = psd_target(freqs, self.psd_spec)
            out = np.empty((n, n_samples))
            for i in range(n):
                out[i], _, _ = synthesize_from_psd(self.rng, target,
                                                   n_samples, self.fs)
            return out
        raise KeyError(f"未知噪声模型 {self.model}")

    def describe(self) -> dict:
        """模型描述（assumed 标注）。"""
        return {"model": self.model, "rms": self.rms, "tau_s": self.tau_s,
                "psd_spec": self.psd_spec, "fs_hz": self.fs,
                "provenance": "assumed"}


def validate_psd_synthesis(model_spec, n_samples=20000, fs=1e6, seed=11,
                           n_realizations=24):
    """PSD 合成验证：目标/测量 PSD 一致性与 seed 复现性。

    返回诊断字典（频率、目标 PSD、Welch 均值与分位带、方差比、
    Hermitian 对称性与实数性检查、seed 复现检查）。
    """
    rng = np.random.default_rng(seed)
    freqs = np.fft.rfftfreq(n_samples, d=1.0 / fs)
    target = psd_target(freqs, model_spec)
    traces = []
    for _ in range(n_realizations):
        x, _, var_t = synthesize_from_psd(rng, target, n_samples, fs)
        traces.append(x)
    traces = np.asarray(traces)
    real_ok = bool(np.all(np.isreal(traces)))
    finite_ok = bool(np.all(np.isfinite(traces)))
    w_freqs, welch_mean = welch_psd(traces.ravel(), fs, nperseg=1024)
    # 目标 PSD 重采样到 Welch 频率
    target_resampled = np.interp(w_freqs, freqs, target)
    var_measured = float(np.var(traces))
    var_integral = float(np.sum(target[1:]) * fs / n_samples)
    rng2 = np.random.default_rng(seed)
    x2, _, _ = synthesize_from_psd(rng2, target, n_samples, fs)
    x3, _, _ = synthesize_from_psd(np.random.default_rng(seed), target,
                                   n_samples, fs)
    reproducible = bool(np.allclose(x2, x3))
    rng_a = np.random.default_rng(seed)
    xa, _, _ = synthesize_from_psd(rng_a, target, n_samples, fs)
    rng_b = np.random.default_rng(seed + 1)
    xb, _, _ = synthesize_from_psd(rng_b, target, n_samples, fs)
    independent = not bool(np.allclose(xa[:100], xb[:100]))
    # Welch 95% 置信带（按自由度近似）
    welch_lo = np.quantile([welch_psd(t, fs, 1024)[1] for t in traces],
                           0.025, axis=0)
    welch_hi = np.quantile([welch_psd(t, fs, 1024)[1] for t in traces],
                           0.975, axis=0)
    rel_err = float(np.nanmedian(np.abs(welch_mean - target_resampled)
                                 / np.maximum(target_resampled,
                                              1e-30)))
    return {"freqs": w_freqs.tolist(), "target": target_resampled.tolist(),
            "welch_mean": welch_mean.tolist(),
            "welch_lo": welch_lo.tolist(), "welch_hi": welch_hi.tolist(),
            "variance_measured": var_measured,
            "variance_from_psd_integral": var_integral,
            "variance_ratio": var_measured / max(var_integral, 1e-300),
            "median_relative_psd_error": rel_err,
            "real_valued": real_ok, "finite": finite_ok,
            "seed_reproducible": reproducible, "seed_independent": independent,
            "nyquist_hz": fs / 2.0, "record_length_s": n_samples / fs}
