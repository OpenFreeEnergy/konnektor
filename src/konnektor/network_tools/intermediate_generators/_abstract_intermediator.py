import abc
import inspect
from collections.abc import Iterator

from gufe import SmallMoleculeComponent
from gufe.tokenization import GufeTokenizable


class Intermediator(GufeTokenizable):
    def __call__(self, *args, **kwargs) -> Iterator:
        return self.generate_intermediate(*args, **kwargs)

    @classmethod
    def _defaults(cls):
        sig = inspect.signature(cls.__init__)

        defaults = {
            param.name: param.default
            for param in sig.parameters.values()
            if param.default is not inspect.Parameter.empty
        }

        return defaults

    @classmethod
    def _from_dict(cls, dct: dict):
        init_args = cls._defaults()
        additional_vas = {}

        for key, value in dct.items():
            if key in init_args:
                init_args.update({key: value})
            else:
                additional_vas[key] = value

        new_obj = cls(**init_args)
        [setattr(new_obj, key, value) for key, value in additional_vas.items()]

        return new_obj

    def _to_dict(self, include_defaults=True) -> dict:
        self_dict = {k: v for k, v in vars(self).items() if not k.startswith("_")}

        if include_defaults:
            add_vars = {k: v for k, v in self.defaults().items() if k not in self_dict}
            self_dict.update(add_vars)

        return self_dict

    @abc.abstractmethod
    def generate_intermediate(
        self, molA: SmallMoleculeComponent, molB: SmallMoleculeComponent
    ) -> Iterator[SmallMoleculeComponent]:
        pass
