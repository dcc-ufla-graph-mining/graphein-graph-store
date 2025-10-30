class IndexedSet:
    def __init__(self, list=None):
        self._list = []
        self._index = {}

        if list:
            for element in list:
                self.add(element)

    def add(self, value):
        if value not in self._index:
            self._index[value] = len(self._list)
            self._list.append(value)
        return self._index[value]

    def contains(self, value):
        return value in self._index
    
    def index(self, value):
        return self._index[value]

    def get(self, idx):
        return self._list[idx] 