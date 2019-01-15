from setuptools import setup

setup(name='two_step_graph_mapper',
      version='0.0.1',
      description='Two Step Graph Mapper',
      url='http://github.com/uio-bmi/two_step_graph_mapper',
      author='Ivar Grytten and Knut Rand',
      author_email='',
      license='MIT',
      zip_safe=False,
      install_requires=['numpy', 'python-coveralls', 'pyvg',
                        'pyfaidx', 'offsetbasedgraph==2.1.3', 'graph_peak_caller', 'tqdm'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points = {
        'console_scripts': ['two_step_graph_mapper=two_step_graph_mapper.command_line_interface:main'],
      }
)