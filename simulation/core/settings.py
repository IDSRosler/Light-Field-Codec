import os
import importlib
from core import global_settings


class ImproperlyConfigured(Exception):
    """Settings is somehow improperly configured"""
    pass


class Settings:
    def __init__(self, settings_module):
        # Load default settings from global_settings
        for setting in dir(global_settings):
            if setting.isupper():
                setattr(self, setting, getattr(global_settings, setting))
        self.SETTINGS_MODULE = settings_module
        mod = importlib.import_module(self.SETTINGS_MODULE)
        for setting in dir(mod):
            if setting.isupper():
                setting_value = getattr(mod, setting)
                setattr(self, setting, setting_value)
        if not self.EXECUTABLE:
            raise ImproperlyConfigured(
                "The EXECUTABLE value must not be empty"
            )
        self.args = []
        self.build_args()

    def build_args(self):
        from itertools import product

        var_args = dict()
        inline_args = []
        ref_args = []
        product_args = []

        for arg, value in self.ARGS.items():
            if type(value) is str:
                if value in self.ARGS:
                    ref_args.append((arg, value))
                else:
                    inline_args.append((arg, value))
            elif type(value) in {int, float, complex}:
                inline_args.append((arg, value))
            elif type(value) in {list, tuple}:
                lst_args = [(arg, x) for x in value]
                var_args.update({arg: lst_args})

        # zip(*args) where args is a list of together values
        for together_args in self.GROUP_TOGETHER_ARGS:
            args = [var_args.pop(a) for a in together_args]
            product_args.append(list(zip(*args)))

        args = []
        product_args.extend(var_args.values())
        for vargs in product(*product_args):
            # Ungroup grouped args.
            flat = []
            for sublist in vargs:
                if type(sublist[0]) in [list, tuple]:
                    flat.extend(sublist)
                else:
                    flat.append(sublist)

            dargs = dict((*inline_args, *flat))
            rargs = [(rarg, dargs[val]) for rarg, val in ref_args]
            sargs = [(sarg, '') for sarg in self.SEQ_ARGS]
            args.append((*dargs.items(), *rargs, *sargs))
        
        for arg in args:
            self.args.append(dict(arg))
        
        

    def parse_variables(self, args):

        to_parse = getattr(self, 'VARIABLES_TO_PARSE', list())
        to_parse.extend(['STDIN', 'STDOUT', 'STDERR'])

        variables = dict(vars(self), **args)

        for var in to_parse:
            value = getattr(self, var)
            if value is not None:
                value = value % variables
                variables[var] = value
        
        return variables
        