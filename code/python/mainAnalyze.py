from StatisticalAnalysis import StatisticalAnalysis

print('started')

#analise direta
obj = StatisticalAnalysis('RESULTADOS-metal-bonds-removed.csv')
obj.analyze()
#obj.analyzePolyhedron('polyhedronValues-version-1.csv') #analise dos poliedros


print('finished')
