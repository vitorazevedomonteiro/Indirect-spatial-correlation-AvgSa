# flake8: noqa
import inspect
from .bakerjayaram08 import corrBJ08
from .duning21 import CrossSpatialCorrDN21
from .jayarambaker09 import SpatialCorrJB09
from .lothbaker13 import CrossSpatialCorrLB13
from .markhvidaEtAl18 import CrossSpatialCorrMCB18
from .monteiroEtAl26 import CrossSpatialCorrMAO26


ALIASES = {
    "bakerjayaram08": corrBJ08,
    "duning21": CrossSpatialCorrDN21,
    "jayarambaker09": SpatialCorrJB09,
    "lothbaker13": CrossSpatialCorrLB13,
    "markhvidaEtAl18": CrossSpatialCorrMCB18,
    "monteiroEtAl26": CrossSpatialCorrMAO26,
}


def select_func_args(function):
    sig = inspect.signature(function)

    func_params = [param.name for param in sig.parameters.values()
                   if param.name != 'self']

    return func_params