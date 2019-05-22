'''Implementation of memoizing decorators
'''
from collections import OrderedDict
from functools import wraps


class MemoizedFunction(object):
    """A memoized function wrapper that manages its own hash table cache.

    The wrapped function's arguments must be hashable.

    Attributes
    ----------
    max_size: int
        The maximum number of distinct invocations to cache
    function: Callable
        The wrapped function
    cache: Mapping
        The hash table cache
    lru: bool
        Whether or not to enforce a LRU eviction strategy

    """
    def __init__(self, function, max_size=100, lru=False):
        self.max_size = max_size
        self.function = function
        self.lru = lru
        self._init_cache()
        wraps(function)(self)

    def _init_cache(self):
        if self.lru:
            self.cache = OrderedDict()
        else:
            self.cache = dict()

    def evict(self):
        if self.lru:
            self.cache.popitem(last=False)
        else:
            self.cache.popitem()

    def __call__(self, *args, **kwargs):
        key = (args, frozenset(kwargs.items()))
        if key not in self.cache:
            if len(self.cache) == self.max_size:
                self.evict()
            self.cache[key] = self.function(*args, **kwargs)
        if self.lru:
            value = self.cache.pop(key)
            self.cache[key] = value
            return value
        return self.cache[key]

    def clear(self):
        self.cache.clear()

    def __len__(self):
        return len(self.cache)

    def __repr__(self):
        template = "{self.__class__.__name__}({self.function}, {size}/{self.max_size})"
        size = len(self.cache)
        return template.format(self=self, size=size)

    @classmethod
    def memoize(cls, maxsize=100):
        """Make a memoization decorator. A negative value of `maxsize` means
        no size limit."""

        def _wrapper(func):
            return cls(func, maxsize)

        return _wrapper


memoize = MemoizedFunction.memoize


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
