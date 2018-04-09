
class SmartDict(dict):
    """Implements Perl's autovivification for dicts.

    Extends dict to create nested dicts if an invalid key is used, rather than throw an exception.

    eg. D = SmartDict()
        print D
        >> {}
        print D[1][2][4]
        >> {1 : {2 : {4 : {} } } }

        E = SmartDict()
        E[1][3] = "Hello World"
        print E[1][3]
        >> Hello World
    """

    def __init__(self):
        dict.__init__(self)
        self.locked = False

    def __missing__(self, key):
        """If given invalid key, store an empty SmartDict object using this key"""
        value = self[key] = type(self)()
        return value
