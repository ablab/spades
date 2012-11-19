
class PathSet:

    def __init__(self, pathset):
        pathset.sort()
        self.pathset = set(tuple(path) for path in pathset)
        self.hash = hash(tuple(self.pathset))

    def __eq__(self, other):
        return self.pathset == other.pathset

    def __hash__(self):
        return self.hash

    def conj(self):
        cpathset = []
        for path in self.pathset:
            cpath = []
            for edge in path[::-1]:
                cpath.append(edge.conj)
            cpathset.append(cpath)
        return PathSet(cpathset)

    def have_common_with(self, other):
        return self.pathset.intersection(other.pathset)

    def crop_left(self):
        cropped = []
        for path in self.pathset:
            cropped.append(path[1:])
        return PathSet(cropped)

    def crop_right(self):
        cropped = []
        for path in self.pathset:
            cropped.append(path[:-1])
        return PathSet(cropped)

    def join_with(self, other):
        pathset = list(self.pathset.union(other.pathset))
        return PathSet(pathset)
