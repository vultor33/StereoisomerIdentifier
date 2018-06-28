from StatisticalAnalysis import StatisticalAnalysis

print('started')

#analise direta
obj = StatisticalAnalysis('Dimetals-Teste.csv')
obj.analyze()
#obj.analyzePolyhedron('Dimetals-Teste.csv'.csv') #analise dos poliedros


print('finished')
