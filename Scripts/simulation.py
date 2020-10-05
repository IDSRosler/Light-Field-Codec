import os

from PlataformaExecucaoMultithread import PlataformaExecucaoMultithread

if __name__ == '__main__':
	lista_comandos = list()

	light_field = {'Bikes', 'Danger_De_Mort', 'Fountain_Vincent2', 'Stone_Pillars_Outside'}
	qps = {0.1, 0.25, 0.5, 1, 2, 5, 7, 10, 15, 20, 50, 100}

	""" Diretório dos arquivos de entrada .ppm """
	input_file = '/media/ubuntu/872f3d86-0680-4dea-b351-493e8c7e7054/home/igor/Git/Full_datasets/'

	""" Diretório de saída """
	output_file = '../output/'

	for name in light_field:
		for qp in qps:
			file_out =  name + '_15_15_13_13_' + str(qp) + '_' + str(qp) + '_' + str(qp) + '_' + str(qp) + '_' + str(qp)

			cmd = dict()
			
			cmd['bin'] = '../cmake-build-debug/./LF_Codec'

			cmd['argv'] = '-input ' + input_file + name + \
						  '/ -blx 15 -bly 15 -blu 13 -blv 13 -qx ' + str(qp) + ' -qy ' + str(qp) + ' -qu ' + str(qp) + ' -qv ' + str(qp) + \
					      ' -lfx 625 -lfy 434 -lfu 13 -lfv 13 -qp ' + str(qp) + \
						  ' -output ' + output_file + file_out + '/'
			
			cmd['output'] = file_out

			lista_comandos.append(cmd)

	plataforma = PlataformaExecucaoMultithread(lista_comandos, output_file)