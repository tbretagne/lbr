from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='lbr',
      version='0.0.1',
      description='Pure Python Implementation of lets_be_rational by Peter Jackel',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Topic :: Office/Business :: Financial',
      ],
      keywords='volatility finance option black',
      url='http://github.com/tbretagne/lbr',
      author='Tanguy Bretagne',
      author_email='tanguy.bretagne@gmail.com',
      license='MIT',
      packages=['lbr'],
      install_requires=[
          'numpy',
      ],
      include_package_data=True,
      zip_safe=False)
