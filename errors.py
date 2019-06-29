class Error(Exception):
    """Base class for other exceptions. User defined exceptions follow"""
    pass

class BadValueError(Error):
    """When the user inputs wrong value"""
    pass


class SimulationTypeError(Error):
    """When the user inputs a wrong simulation type"""
    pass

class StarsmasherInternalError(Error):
    """When Starsmasher reports some internal error. Check Log"""
    pass


