from setuptools import setup

requires = [
    'pyramid',
    'waitress',
]

setup(
    name='ui',
    install_requires=requires,
    entry_points={
        'paste.app_factory': [
            'main = ui:main'
        ],
    },
)
