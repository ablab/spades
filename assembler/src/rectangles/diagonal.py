import math
import utils
import experimental

class Key(object):
    def __init__(self, first, second, pathset):
        self.first = first
        self.second = second
        self.pathset = pathset
        self.hash = hash((first, second))

    def __eq__(self, other):
        if self.first != other.first or self.second != other.second:
            return False
        if experimental.filter == experimental.Filter.pathsets:
            if not self.pathset.have_common_with(other.pathset):
                return False
        return True

    def __hash__(self):
        return self.hash

    def join_with(self, other):
        assert self.first == other.first
        assert self.second == other.second
        if experimental.filter == experimental.Filter.pathsets:
            pathset = self.pathset.join_with(other.pathset)
        else:
            pathset = None
        return Key(self.first, self.second, pathset)


class Diagonal(object):

    def __init__(self, rectangle, d, D, pathset):
        self.rectangle = rectangle
        self.D = D
        self.prd_support = 0.0
        self.offseta, self.offsetb, self.offsetc, self.offsetd = Diagonal.offsets(rectangle.e1, rectangle.e2, d, D)
        assert self.offseta != self.offsetc
        assert self.offsetb != self.offsetd
        self.pathset = pathset
        self.d = d
        self.key1 = Diagonal.key(self.rectangle.e1, self.offseta, self.rectangle.e2, self.offsetb, pathset)
        self.key2 = Diagonal.key(self.rectangle.e1, self.offsetc, self.rectangle.e2, self.offsetd, pathset)

    def inc_closest(self, D, weight):
        r = self.rectangle
        # TODO: optimize (on top level)
        if abs(self.D - D) == abs(r.get_closest_diagonal(D).D - D):
            # TODO: don't count point between diagonals twice
            self.prd_support += weight

    def inc_in_range(self, D, weight, delta):
        if abs(self.D - D) <= delta:
            self.prd_support += weight

    def inc_normal(self, D, weight, var):
        # probability being normal
        A = abs(self.D - D) - 0.5
        B = A + 1
        probability = 0.5 * (utils.erf(B / var / math.sqrt(2)) - utils.erf(A / var / math.sqrt(2)))
        self.prd_support += probability * weight

    def inc_hist(self, D, weight, hist):
        x = self.d + int(round(D - self.D)) # TODO: check sign!
        if x in hist:
            self.prd_support += hist[x] * weight

    def inc_prd_support(self, D, weight, delta, config):
        if weight == 0:
            assert self.rectangle.e1 == self.rectangle.e2
            return
        if experimental.filter == experimental.Filter.spades:
            self.inc_in_range(D, weight, delta)
        elif experimental.filter == experimental.Filter.closest:
            self.inc_closest(D, weight)
        elif experimental.filter == experimental.Filter.normal:
            self.inc_normal(D, weight, config.is_var)
        elif experimental.filter == experimental.Filter.hist:
            self.inc_hist(D, weight, config.hist)
        else:
            assert False

    def __repr__(self):
        return 'D(%s|%s,%d)' % (str(self.rectangle.e1), str(self.rectangle.e2), self.D)

    def support(self):
        if self.D == 0:
            return 1e100
        else:
            return self.prd_support / (self.offsetc - self.offseta)

    @staticmethod
    def key(e1, offset1, e2, offset2, pathset):
        # key = ((e, o), v, pathset) or (v, (e, o), pathset) or (v, v, pathset)
        first = e1.v1 if offset1 == 0 else ( e1.v2 if offset1 == e1.len else (e1, offset1) )
        second = e2.v1 if offset2 == 0 else ( e2.v2 if offset2 == e2.len else (e2, offset2) )
        if experimental.filter == experimental.Filter.pathsets:
            if offset2 == 0:
                pathset = pathset.crop_right()
            if offset1 == e1.len:
                pathset = pathset.crop_left()
        return Key(first, second, pathset)

    @staticmethod
    def offsets(e1, e2, d, D):
        l1, l2 = e1.len, e2.len
        # calculate offsets in rectangle
        if d >= D:
            offseta, offsetb = 0, d - D
        else:
            offseta, offsetb = D - d, 0
        if d >= D + l2 - l1:
            offsetc, offsetd = l2 + (D - d), l2
        else:
            offsetc, offsetd = l1, l1 + (d - D)
        assert offsetc - offseta == offsetd - offsetb, "Should be main diagonal (it's a bug)"
        assert 0 <= offseta <= l1, (offseta, l1)
        assert 0 <= offsetb <= l2, (offsetb, l2)
        assert 0 <= offsetc <= l1, (offsetc, l1)
        assert 0 <= offsetd <= l2, (offsetd, l2)
        return offseta, offsetb, offsetc, offsetd
