try:
    from importlib.metadata import version, PackageNotFoundError
    try:
        __version__ = version("xcid")
    except PackageNotFoundError:
        from ._version import __version__  # type: ignore
except ImportError:
    __version__ = None

from .xcid import XCID

__all__ = ["XCID", "__version__"]
