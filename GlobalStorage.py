from ACTORS_SETUP import *

@ray.remote
class GlobalStorage:
    def __init__(self):
        self.store = {}

    def put(self, key, value):
        self.store[key] = value

    def get(self, key):
        return self.store.get(key, None)

    def return_all(self):
        return self.store