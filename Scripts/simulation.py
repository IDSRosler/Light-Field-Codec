import os

from PlataformaExecucaoMultithread import PlataformaExecucaoMultithread

if __name__ == '__main__':
	lista_comandos = list()

	light_field = ['Bikes']
	QUANTIZATIONS = [3]
	TRANSFORMS = ['DCT', 'DST_II', 'DST_VII']
	""" Diretório dos arquivos de entrada .ppm """
	input_file = '/home/cristian/ufpel/lfcodec/datasets/'

	""" Diretório de saída """
	output_file = '/home/cristian/ufpel/lfcodec/results_luma/'

	for name in light_field:
		for qp in QUANTIZATIONS:
			for tx in TRANSFORMS:
				file_out =  name + '_15_15_13_13_' + str(qp) + '_' + str(qp) + '_' + str(qp) + '_' + str(qp) + '_' + str(qp) + '_' + tx

				cmd = dict()
				
				cmd['bin'] = '../build/LF_Codec'

				cmd['argv'] = '-input ' + input_file + name + \
							'/ -blx 15 -bly 15 -blu 13 -blv 13 -qx ' + str(qp) + ' -qy ' + str(qp) + ' -qu ' + str(qp) + ' -qv ' + str(qp) + \
							' -lfx 625 -lfy 434 -lfu 13 -lfv 13 -qp ' + str(qp) + ' -transform ' + tx + \
							' -output ' + output_file + file_out + '/' 
				
				cmd['output'] = file_out

				lista_comandos.append(cmd)

	plataforma = PlataformaExecucaoMultithread(lista_comandos, output_file)