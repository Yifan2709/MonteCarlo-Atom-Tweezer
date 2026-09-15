"""有解析预期的多站点 IRB 回归，防止跨序列状态及生存标签串写。"""
from types import SimpleNamespace

import numpy as np
import pytest

from level5_finite_pulse_irb.interleaved_rb import run_interleaved_rb
from level5_finite_pulse_irb.rb_fitting import irb_estimate


@pytest.mark.parametrize("n_sites", [1, 3, 47])
def test_move_updates_every_site_of_only_selected_sequence(n_sites):
    starts = []

    class IdentityGates:
        def apply(self, psi, cliff, trace, tstart, *args, **kwargs):
            starts.append(tstart.copy())

        def clifford_duration_s(self, index):
            return 0.25

    class BitFlipMove:
        def propagate(self, psi, *args, **kwargs):
            psi[:] = psi[:, ::-1]

    counts = [0, 1, 2, 60]
    seqs = [{"M": m, "cliffords": [0] * 60, "inverse": 0,
             "string_id": 0, "sequence_index": 0} for m in counts]
    pool = SimpleNamespace(n_traj=1, roundtrip_success=np.array([False]),
                           total_duration_s=2.0)
    result = run_interleaved_rb(seqs, BitFlipMove(), IdentityGates(),
                                pool, pool, 0.0, 0.0, seed=19,
                                n_sites=n_sites)
    expected = np.repeat([[1.0], [0.0], [1.0], [1.0]], n_sites, axis=1)
    np.testing.assert_allclose(result["conditional_return"], expected)
    np.testing.assert_array_equal(result["survival"],
                                  np.repeat([[1], [0], [0], [0]], n_sites, axis=1))
    # 最后的 inverse 必须看到实际流逝的 gate + move 时间。
    np.testing.assert_allclose(starts[-1].reshape(4, n_sites),
                               np.repeat((15 + 2 * np.array(counts))[:, None],
                                         n_sites, axis=1))


def test_fixed_n_m_scan_cannot_use_standard_irb_ratio():
    est = irb_estimate(0.9967, 0.995, design="fixed_cliffords_variable_moves")
    assert est["interpretable"] is False
    assert "f_avg_interleaved" not in est


@pytest.mark.parametrize("p_ref,p_int", [(0.9967, 1.0), (np.nan, 0.9),
                                        (0.9, np.inf), (1.1, 0.9)])
def test_invalid_irb_fidelity_is_not_published(p_ref, p_int):
    est = irb_estimate(p_ref, p_int)
    assert est["interpretable"] is False
    assert "f_avg_interleaved" not in est
