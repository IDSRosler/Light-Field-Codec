import os

from PlataformaExecucaoMultithread import PlataformaExecucaoMultithread

if __name__ == '__main__':
	lista_comandos = list()

	"""light_field = {'Bikes', 'DangerDeMort', 'FountainVincent2', 'StonePillarsOutside'}
	qps = {0.1, 0.25, 0.5, 1, 2, 5, 7, 10, 15, 20, 50, 100}"""

	light_field = {'Bikes', 'DangerDeMort', 'FountainVincent2', 'StonePillarsOutside'}
	qps = {0.1, 1, 5, 10}

	""" Diretório dos arquivos de entrada .ppm """
	input_file = '/home/igor/Documentos/Git/Full_datasets/'

	""" Diretório de saída """
	output_file = '/home/igor/Documentos/Git/Outputs/output_Last_Run_0_1/'

	""" Tamanho da partição de bloco """
	block_particion = '4x4x4x4'

	for name in light_field:
		for qp in qps:
			#file_out =  name + '_15_15_13_13_' + str(qp) + '_' + str(qp) + '_' + str(qp) + '_' + str(qp) + '_' + str(qp)
			file_out =  name + '_16_16_16_16_' + str(qp) + "_13_13"

			cmd = dict()
			
			cmd['bin'] = '../cmake-build-debug/./Light_Field_Codec'

			cmd['argv'] = '-input ' + input_file + name + \
						  '/ -blx 15 -bly 15 -blu 13 -blv 13 -qx ' + str(qp) + ' -qy ' + str(qp) + ' -qu ' + str(qp) + ' -qv ' + str(qp) + \
					      ' -lfx 625 -lfy 434 -lfu 13 -lfv 13 -qp ' + str(qp) + \
						  ' -output ' + output_file + file_out + '/'
			
			cmd['output'] = file_out

			lista_comandos.append(cmd)

	plataforma = PlataformaExecucaoMultithread(lista_comandos, output_file)