try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup


with open('README.md') as f:
    long_description = f.read()
    
setup(
    name='k2gap',
    version=1.0.0,
    description='A module for the K2GAP target selection function',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/sanjibs/k2gap',
    author='Sanjib Sharma',
    author_email='bug.sanjib@gmail.com',
    license='New BSD',
    classifiers=[
        'Development Status :: 6 - Mature',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3.7',
    ],
    install_requires=['numpy','json'],
    packages=['k2gap'],
    package_data={'': ['k2circles.json','AUTHORS.rst','README.md','LICENSE']},
)
