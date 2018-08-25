from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import setuptools


class GetPybindInclude(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ext_modules = [
    Extension(
        'nadavca.dtw',
        sources=['nadavca/dtw/dtwmodule.cpp',
                 'nadavca/dtw/kmer_model.cpp',
                 'nadavca/dtw/probability.cpp',
                 'nadavca/dtw/sequence.cpp',
                 'nadavca/dtw/node.cpp',
                 'nadavca/dtw/dtw.cpp'],
        include_dirs=[
            'nadavca/dtw/',
            # Path to pybind11 headers
            GetPybindInclude(),
            GetPybindInclude(user=True)
        ],
        language='c++'
    ),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as file:
        file.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([file.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    def build_extensions(self):
        compiler_type = self.compiler.compiler_type
        opts = self.c_opts.get(compiler_type, [])
        if compiler_type == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append('-std=c++11')
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif compiler_type == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)


setup(
    name='nadavca',
    version='0.1',
    author='Eduard Batmendijn',
    author_email='ebatmendijn@gmail.com',
    url='https://github.com/baklazan/nadavca',
    description='NAnopore DAta Variant CAller',
    install_requires=['pybind11>=2.2', 'simplesam', 'h5py', 'numpy'],
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},
    scripts=['bin/nadavca'],
    packages=['nadavca'],
    include_package_data=True,
    zip_safe=False
)
