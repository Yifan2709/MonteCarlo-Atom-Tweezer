"""Typed configuration errors with stable codes and field paths.

Every user-facing validation failure raises ToolkitError carrying a stable
error code (E_*) and the dotted field path of the offending input, so callers
can locate the problem without parsing prose. Codes are part of the public
contract: new codes may be added, existing meanings never silently change.
"""
from __future__ import annotations


class ToolkitError(ValueError):
    """Configuration or contract error with a stable code and field path."""

    def __init__(self, code: str, message: str, field: str = ""):
        super().__init__(f"[{code}] {f'{field}: ' if field else ''}{message}")
        self.code = code
        self.field = field
        self.message = message


# Stable error codes (public contract).
E_UNKNOWN_FIELD = "E_UNKNOWN_FIELD"          # key present but not in schema
E_MISSING_FIELD = "E_MISSING_FIELD"          # required key absent
E_BAD_TYPE = "E_BAD_TYPE"                    # value has wrong type
E_NOT_FINITE = "E_NOT_FINITE"                # NaN/inf where finite required
E_BAD_VALUE = "E_BAD_VALUE"                  # finite but out of domain
E_UNIT = "E_UNIT"                            # unparseable/unsupported unit
E_TIME_ORDER = "E_TIME_ORDER"                # duplicate/decreasing times
E_LENGTH_MISMATCH = "E_LENGTH_MISMATCH"      # array length disagreement
E_MODEL_UNKNOWN = "E_MODEL_UNKNOWN"          # named model not implemented
E_CAL_RANGE = "E_CAL_RANGE"                  # calibration range violated
E_NOT_DIVISIBLE = "E_NOT_DIVISIBLE"          # grid divisibility requirement
E_UNSUPPORTED = "E_UNSUPPORTED"              # valid request outside supported scope
E_CONFLICT = "E_CONFLICT"                    # fields inconsistent with each other
