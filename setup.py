from setuptools import setup, find_packages

import versioneer

setup(
    name="q2-gamma",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={
        'q2_gamma.plot': ['assets/dist/*', 'assets/dist/licenses/*']
    },
    author="Evan Bolyen",
    author_email="ebolyen@gmail.com",
    description="TBD",
    license='BSD-3-Clause',
    url="TBD",
    entry_points={
        'qiime2.plugins': ['q2-gamma=q2_gamma.plugin_setup:plugin']
    },
    zip_safe=False,
)

