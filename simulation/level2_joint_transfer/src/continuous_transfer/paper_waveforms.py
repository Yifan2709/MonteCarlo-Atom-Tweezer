"""Fig. 6b constant-jerk cubic; historical S-curve remains selectable."""
import numpy as np
from level2_joint_transfer.waveforms import PickupSequentialWaveform, _out


class PaperManualWaveform(PickupSequentialWaveform):
    """x/d=3s²-2s³: zero endpoint velocity, constant interior jerk.

    Acceleration jumps at the start/end of motion. It must not be silently
    replaced by an S-curve with zero endpoint acceleration.
    """
    def _terms(self, t):
        u = self._move_u(t)
        span = self.final_center_m if self.direction == "pickup" else -self.final_center_m
        offset = 0.0 if self.direction == "pickup" else self.final_center_m
        tau = self.move_duration_s
        active = (u > 0) & (u < 1)
        return (offset+span*(3*u*u-2*u*u*u), span*6*u*(1-u)/tau,
                np.where(active,span*(6-12*u)/tau**2,0),
                np.where(active,-12*span/tau**3,0))

    def center(self,t): return _out(self._terms(t)[0],t)
    def center_velocity(self,t): return _out(self._terms(t)[1],t)
    def center_acceleration(self,t): return _out(self._terms(t)[2],t)
    def center_jerk(self,t): return _out(self._terms(t)[3],t)

    def describe(self):
        return dict(super().describe(), position_profile="paper_single_cubic_constant_jerk")
