# misc untils

from math import log

class ColumnDensity:
    def __init__(self, start):
        self.data = float(start)
    def __str__(self):
        return str(self.data)
    def __add__(self, other):
        return ColumnDensity(log(10**self.data+10**other.data, 10))
    def __sub__(self, other):
        return ColumnDensity(log(10**self.data-10**other.data, 10))
    def __mul__(self, other):
        return ColumnDensity(self.data+other.data)
    def __div__(self, other):
        return ColumnDensity(self.data-other.data)
    def __cmp__(self, other):
        if self.data < other.data:
            return -1
        elif self.data == other.data:
            return 0
        else:
            return 1
