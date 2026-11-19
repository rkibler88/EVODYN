"""Exception hierarchy.

EVODYN distinguishes three failure modes, because they mean different things to
a user: a configuration error is the user's to fix, a missing external tool is
an installation problem, and a mapping error is a correctness failure that must
never be silently tolerated.
"""


class EvodynError(Exception):
    """Base class for all EVODYN errors."""


class ConfigError(EvodynError):
    """The run configuration is incomplete or self-inconsistent."""


class ExternalToolError(EvodynError):
    """A required external executable is absent or exited non-zero."""


class MappingError(EvodynError):
    """A residue-numbering or coordinate mapping assertion failed.

    This is deliberately fatal. The EGFR study found numbering errors that had
    propagated into figures and text (residue 853 labelled as the DFG motif,
    which is actually D855-F856-G857); the corresponding checks here stop the
    run rather than annotate the output.
    """


class StageSkipped(EvodynError):
    """A stage could not run and was skipped; recorded in the run manifest."""
