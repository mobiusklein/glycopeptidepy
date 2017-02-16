from collections import Iterable

sentinel = object()


class peekable(Iterable):
    sentinel = sentinel

    def __init__(self, iterable):
        self.iterable = iter(iterable)
        try:
            self._lookahead = next(self.iterable)
        except StopIteration:
            self._lookahead = sentinel

    def __iter__(self):
        # Counter to determine if we've passed through the
        # inner loop at all
        i = 0
        # Default stopping value for x to make accidental reuse
        # of the iterator safer
        x = sentinel
        for x in self.iterable:
            val = self._lookahead
            self._lookahead = x
            yield (val, x)
            i += 1
            x = sentinel
        else:
            val = self._lookahead
            self._lookahead = x
        if i > 0:
            yield val, x
        return

    @property
    def peek(self):
        return self._lookahead
