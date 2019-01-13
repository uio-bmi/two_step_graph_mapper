from setuptools import setup

setup(name='rough_graph_mapper',
      version='0.0.1',
      description='Rough Graph Mapper',
      url='http://github.com/uio-bmi/rough_graph_mapper',
      author='Ivar Grytten and Knut Rand',
      author_email='',
      license='MIT',
      zip_safe=False,
      install_requires=['numpy', 'python-coveralls',
                        'pyfaidx', 'offsetbasedgraph', 'graph_peak_caller', 'tqdm'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points = {
        'console_scripts': ['rough_graph_mapper=rough_graph_mapper.command_line_interface:main'],
      }

      )