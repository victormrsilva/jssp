instances = ['vi_33_28', 'vi_44_39', 'ft06']
horizon = [33, 44, 60]

for i in range(len(instances)):
    for select in range(4):
        for parameter in range(3, 9):
            for gomory in range(2):
                for zero in range(2):
                    for mir in range(2):
                        l = parameter
                        if select == 1:
                            l += 2
                        str = '{}{}{}{}{}'.format(select, parameter, gomory, zero, mir)
                        filename = '{}_{}.cfg'.format(instances[i], str)
                        with open(filename, 'w') as f:
                            f.write('instance_name = {}\n'.format(instances[i]))
                            f.write('horizon = {}\n'.format(horizon[i]))
                            f.write('problem = 1\n')
                            f.write('mip_general_cliques_select = {}\n'.format(select))
                            f.write('mip_general_cliques_parameter = {}\n'.format(l))
                            f.write('cut_type_gomory = {}\n'.format(gomory))
                            f.write('cut_type_gomory_limit = 30\n')
                            f.write('cut_type_zero_half = {}\n'.format(zero))
                            f.write('cut_type_mir = {}\n'.format(mir))
                            f.write('cut_mip_maximum = 1000\n')
                        f.close()
                        with open('teste.sh', 'a') as f:
                            f.write('python3 jssp_compact.py testes/{} > testes/{}.log\n'.format(filename, filename))
                        f.close()
