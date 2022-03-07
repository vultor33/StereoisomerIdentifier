from collections import defaultdict
from identifier.ErrorMessages import ErrorMessages

class Graph(object):
	""" Graph data structure, undirected by default. """
	def __init__(self, connections, directed=False):
		self.__errorMessages_ = ErrorMessages()
		self._graph = defaultdict(set)
		self._directed = directed
		self.add_connections(connections)
		
		self.__maxGraphTries = 100000
		self.__graphFindTries = 0

	def add_connections(self, connections):
		""" Add connections (list of tuple pairs) to graph """
		for node1, node2 in connections:
			self.add(node1, node2)

	def add(self, node1, node2):
		""" Add connection between node1 and node2 """

		self._graph[node1].add(node2)
		if not self._directed:
			self._graph[node2].add(node1)

	def remove(self, node):
		""" Remove all references to node """
		if node in self._graph:
			del self._graph[node]
		
		for key in self._graph:
			if node in self._graph[key]:
				self._graph[key].remove(node)

	def getNodeConnections(self, node):
		return list(self._graph[node])

	def isConnected(self, node1, node2):
		""" Is node1 connected to node2 """
		self.__graphFindTries = 0
		anyPath = self._findPath(node1, node2)
		return not anyPath is None
		
	def _findPath(self, node1, node2, path=[]):
		""" Find any path between node1 and node2 (may not be shortest) """

		self.__graphFindTries +=1
		if self.__graphFindTries > self.__maxGraphTries:
			raise Exception(self.__errorMessages_.getGraphError())

		path = path + [node1]
		if node1 == node2:
			return path
		if node1 not in self._graph:
			return None
		for node in self._graph[node1]:
			if node not in path:
				new_path = self._findPath(node, node2, path)
				if new_path:
					return new_path
		
		return None

	def __str__(self):
		return '{}({})'.format(self.__class__.__name__, dict(self._graph))
