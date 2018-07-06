from StatisticalAnalysis import StatisticalAnalysis

print('started')

#analise direta
obj = StatisticalAnalysis('RESULTADOS-04-07.csv')
obj.analyze()
#obj.analyzePolyhedron('Dimetals-Teste.csv'.csv') #analise dos poliedros


print('finished')
