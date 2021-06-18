import sys
import os
import re
import traceback

from setuptools import setup, find_packages, Extension

from distutils.command.build_ext import build_ext
from distutils.errors import (CCompilerError, DistutilsExecError,
                              DistutilsPlatformError)


def has_option(name):
    try:
        sys.argv.remove('--%s' % name)
        return True
    except ValueError:
        pass
    # allow passing all cmd line options also as environment variables
    env_val = os.getenv(name.upper().replace('-', '_'), 'false').lower()
    if env_val == "true":
        return True
    return False


include_diagnostics = has_option("include-diagnostics")
force_cythonize = has_option("force-cythonize")


def make_extensions():
    is_ci = bool(os.getenv("CI", ""))
    macros = []
    try:
        from Cython.Build import cythonize
        cython_directives = {
            'embedsignature': True,
            "profile": include_diagnostics
        }
        if include_diagnostics:
            macros.append(("CYTHON_TRACE_NOGIL", "1"))
        if is_ci and include_diagnostics:
            cython_directives['linetrace'] = True
        extensions = cythonize([
            Extension("glycopeptidepy._c.collectiontools", sources=['glycopeptidepy/_c/collectiontools.pyx']),
            Extension("glycopeptidepy._c.count_table", sources=['glycopeptidepy/_c/count_table.pyx']),
            Extension("glycopeptidepy._c.parser", sources=['glycopeptidepy/_c/parser.pyx']),
            Extension("glycopeptidepy._c.structure.base", sources=['glycopeptidepy/_c/structure/base.pyx']),
            Extension("glycopeptidepy._c.structure.fragment", sources=['glycopeptidepy/_c/structure/fragment.pyx']),
            Extension("glycopeptidepy._c.structure.constants", sources=['glycopeptidepy/_c/structure/constants.pyx']),
            Extension("glycopeptidepy._c.structure.glycan", sources=[
                      'glycopeptidepy/_c/structure/glycan.pyx']),
            Extension("glycopeptidepy._c.structure.sequence_methods",
                      sources=['glycopeptidepy/_c/structure/sequence_methods.pyx']),
            Extension("glycopeptidepy._c.structure.modification.rule",
                      sources=['glycopeptidepy/_c/structure/modification/rule.pyx']),
            Extension("glycopeptidepy._c.structure.modification.modification",
                      sources=['glycopeptidepy/_c/structure/modification/modification.pyx']),
            Extension("glycopeptidepy._c.structure.modification.source",
                      sources=['glycopeptidepy/_c/structure/modification/source.pyx']),
            Extension("glycopeptidepy._c.structure.fragmentation_strategy.base",
                      sources=["glycopeptidepy/_c/structure/fragmentation_strategy/base.pyx"]),
            Extension("glycopeptidepy._c.structure.fragmentation_strategy.peptide",
                      sources=["glycopeptidepy/_c/structure/fragmentation_strategy/peptide.pyx"]),
            Extension("glycopeptidepy._c.structure.fragmentation_strategy.glycan",
                      sources=["glycopeptidepy/_c/structure/fragmentation_strategy/glycan.pyx"]),
            Extension("glycopeptidepy._c.algorithm", sources=[
                      'glycopeptidepy/_c/algorithm.pyx']),
        ], compiler_directives=cython_directives, force=force_cythonize)
    except ImportError:
        extensions = [
            Extension("glycopeptidepy._c.collectiontools", sources=['glycopeptidepy/_c/collectiontools.c']),
            Extension("glycopeptidepy._c.count_table", sources=['glycopeptidepy/_c/count_table.c']),
            Extension("glycopeptidepy._c.parser", sources=['glycopeptidepy/_c/parser.c']),
            Extension("glycopeptidepy._c.structure.base", sources=['glycopeptidepy/_c/structure/base.c']),
            Extension("glycopeptidepy._c.structure.fragment", sources=['glycopeptidepy/_c/structure/fragment.c']),
            Extension("glycopeptidepy._c.structure.constants", sources=['glycopeptidepy/_c/structure/constants.c']),
            Extension("glycopeptidepy._c.structure.glycan", sources=[
                      'glycopeptidepy/_c/structure/glycan.c']),
            Extension("glycopeptidepy._c.structure.sequence_methods",
                      sources=['glycopeptidepy/_c/structure/sequence_methods.c']),
            Extension("glycopeptidepy._c.structure.modification.rule",
                      sources=['glycopeptidepy/_c/structure/modification/rule.c']),
            Extension("glycopeptidepy._c.structure.modification.modification",
                      sources=['glycopeptidepy/_c/structure/modification/modification.c']),
            Extension("glycopeptidepy._c.structure.modification.source",
                      sources=['glycopeptidepy/_c/structure/modification/source.c']),
            Extension("glycopeptidepy._c.structure.fragmentation_strategy.base",
                      sources=["glycopeptidepy/_c/structure/fragmentation_strategy/base.c"]),
            Extension("glycopeptidepy._c.structure.fragmentation_strategy.peptide",
                      sources=["glycopeptidepy/_c/structure/fragmentation_strategy/peptide.c"]),
            Extension("glycopeptidepy._c.structure.fragmentation_strategy.glycan",
                      sources=["glycopeptidepy/_c/structure/fragmentation_strategy/glycan.c"]),
            Extension("glycopeptidepy._c.algorithm", sources=['glycopeptidepy/_c/algorithm.c']),
        ]
    return extensions


ext_errors = (CCompilerError, DistutilsExecError, DistutilsPlatformError)
if sys.platform == 'win32':
    # 2.6's distutils.msvc9compiler can raise an IOError when failing to
    # find the compiler
    ext_errors += (IOError,)


class BuildFailed(Exception):

    def __init__(self):
        self.cause = sys.exc_info()[1]  # work around py 2/3 different syntax

    def __str__(self):
        return str(self.cause)


class ve_build_ext(build_ext):
    # This class allows C extension building to fail.

    def run(self):
        try:
            build_ext.run(self)
        except DistutilsPlatformError:
            traceback.print_exc()
            raise BuildFailed()

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except ext_errors:
            traceback.print_exc()
            raise BuildFailed()
        except ValueError:
            # this can happen on Windows 64 bit, see Python issue 7511
            traceback.print_exc()
            if "'path'" in str(sys.exc_info()[1]):  # works with both py 2/3
                raise BuildFailed()
            raise


cmdclass = {}

cmdclass['build_ext'] = ve_build_ext


def status_msgs(*msgs):
    print('*' * 75)
    for msg in msgs:
        print(msg)
    print('*' * 75)


def version():
    return re.sub(r'[\s\'"\n]', '', open("glycopeptidepy/version.py").readline().split("=")[1])


required = []
with open('requirements.txt') as f:
    required = f.read().splitlines()


def run_setup(include_cext=True):
    setup(
        name='glycopeptidepy',
        version=version(),
        packages=find_packages(),
        install_requires=required,
        include_package_data=True,
        ext_modules=make_extensions() if include_cext else None,
        cmdclass=cmdclass,
        zip_safe=False,
        package_data={
            'glycopeptidepy': ["*.csv", "*.xml", "*.json", "data/*.csv"],
            'glycopeptidepy.structure': [
                "structure/modification/data/*.csv",
                "structure/modification/data/*.json"]
        },
    )


try:
    run_setup(True)
except Exception as exc:
    run_setup(False)

    status_msgs(
        "WARNING: The C extension could not be compiled, " +
        "speedups are not enabled.",
        "Plain-Python build succeeded."
    )
