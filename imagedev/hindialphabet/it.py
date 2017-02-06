import numpy as np
class iv(object):
    def merge(self,iv1):
        a = [self[0],self[1],iv1[0],iv1[1]]
        a = sorted(a)
        return [[a[0],a[1]],[a[2],a[3]]]

class it(object):

    def __init__(lists):
        for list in lists:
            self.append(list)
