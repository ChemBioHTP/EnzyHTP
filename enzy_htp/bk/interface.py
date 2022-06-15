from .molecular_mechanics import AmberInterface


class Interface:
    def __init__(self, config):
        self._amber = AmberInterface(config)

    def amber(self):
        return self._amber

    pass
