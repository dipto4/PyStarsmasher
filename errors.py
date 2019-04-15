class Error(Exception):
    """Base class for other exceptions. User defined exceptions follow"""
    pass

class ValueTooSmallError(Error):
    """When the user inputs values smaller than allowed value"""
    pass

class ValueTooLargeError(Error):
    """When the user inputs values larger than allowed value"""
    pass

class SimulationTypeError(Error):
    """When the user inputs a wrong simulation type"""
    pass

class StarsmasherInternalError(Error):
    """When Starsmasher reports some internal error. Check Log"""
    pass


