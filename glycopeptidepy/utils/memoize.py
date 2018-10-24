from functools import wraps


def memoize(maxsize=100):
    """Make a memoization decorator. A negative value of `maxsize` means
    no size limit."""
    def deco(f):
        """Memoization decorator. Items of `kwargs` must be hashable."""
        memo = {}

        @wraps(f)
        def func(*args, **kwargs):
            key = (args, frozenset(kwargs.items()))
            if key not in memo:
                if len(memo) == maxsize:
                    memo.popitem()
                memo[key] = f(*args, **kwargs)
            return memo[key]
        func.memo = memo
        return func
    return deco


class FragmentCachingMixin(object):
    '''Implements a cache over get_fragments() that saves the last arguments and result
    value, and returns it instead of recomputing result if the next invocation's arguments
    match.
    '''
    def __init__(self, *args, **kwargs):
        self.__fragments_value = None
        self.__fragments_arguments = None
        super(FragmentCachingMixin, self).__init__(*args, **kwargs)

    def get_fragments(self, *args, **kwargs):
        if self.__fragments_arguments == (args, kwargs):
            return self.__fragments_value
        else:
            v = super(FragmentCachingMixin, self).get_fragments(*args, **kwargs)
            self.__fragments_value = tuple(v)
            self.__fragments_arguments = (args, kwargs)
            return self.__fragments_value

    def clear_fragments_cache(self):
        self.__fragments_value = None
        self.__fragments_arguments = None

    @classmethod
    def implement(cls, tp):
        return type(tp.__name__, (cls, tp), {})
