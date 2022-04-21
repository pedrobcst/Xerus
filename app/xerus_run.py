from typing import List, Union

from Xerus import XRay


def run_xerus(args_xerus: dict, args_analysis: dict) -> XRay:
    xerus_object = XRay(**args_xerus)
    xerus_object.analyze(**args_analysis)
    return xerus_object

def run_opt(xerus_object: XRay, index_list: Union[int, List[int]], args_opt: dict) -> XRay:
    r = xerus_object
    r.initialize_optimizer(index_list)
    r.run_optimizer(**args_opt)
    return r
