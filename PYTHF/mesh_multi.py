# AUTHOR: CJC
# CREATE TIME: 2022/11/3
# DESCRIPTION: SIMPLE SCRIPT


from mesh_prop import *
from bisect import bisect_left


class MeshMulti:
    def __init__(self,
                 filename,
                 ini_P,
                 ini_T,
                 porosity,
                 Shyd,
                 Sgas,
                 xlen=[],
                 k=[0.0, 0.0, 0.0],
                 ztop=[-0.5, -0.5, -0.5],
                 prop_p=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                 prop_h=[1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0],
                 prop_g=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
                 chemical=0.0,
                 dict_material=None,
                 dict_eos=None,
                 ):
        self.file = filename
        self.ini_P = ini_P
        self.ini_T = ini_T
        self.poro1 = porosity
        self.Shyd = Shyd
        self.Sgas = Sgas
        self.k = k
        self.ztop = ztop
        self.poro2 = prop_p
        self.h = prop_h
        self.g = prop_g
        self.c = chemical
        self.zmax = max(ztop)
        self.xlen = xlen
        self.xposition = np.array(xlen).cumsum()
        self.material = dict_material if dict_material else material
        self.eos = dict_eos if dict_eos else eos

    @staticmethod
    def location(line):

        l_x = float(line[50:60])
        l_y1 = line[60:70]
        if l_y1 == '          ':
            l_y = 0
        else:
            l_y = l_y1
        l_z = float(line[70:80])
        return l_x, l_y, l_z

    @staticmethod
    def is_number(s):
        """
        judge if string is number
        :param s: selected str part
        :return: boolean value
        """
        try:
            float(s)
            return True
        except ValueError:
            pass

        try:
            import unicodedata
            unicodedata.numeric(s)
            return True
        except (TypeError, ValueError):
            pass

        return False

    def set_material(self, line, rocks_id, active=True):
        """

        :param active: if element is active
        :param line: str of a single line
        :param rocks_id: material label
        :return: modified line
        """

        if self.is_number(line[15:20]):

            if rocks_id != 3:

                if active:

                    character = list(line)
                    character[15:20] = '    ' + str(rocks_id)
                    c2 = ''.join(character)
                    return c2
                else:

                    character = list(line.strip())
                    character[15:20] = '    ' + str(rocks_id)
                    c2 = ''.join(character)
                    c3 = c2 + ' I' + '\n'
                    return c3

            else:

                character = list(line.strip())
                character[15:20] = '    ' + str(rocks_id)
                c2 = ''.join(character)
                c3 = c2 + ' I' + '\n'
                return c3

        else:

            if rocks_id != 3:

                character = list(line)
                label = [k for k, v in self.material.items() if v == rocks_id]
                character[15:20] = label
                c2 = ''.join(character)
                return c2

            else:

                character = list(line.strip())
                label = [k for k, v in self.material.items() if v == rocks_id]
                character[15:20] = label
                c2 = ''.join(character)
                c3 = c2 + ' I' + '\n'
                return c3

    def set_subdomain(self,
                      filename2,
                      bound=[-100000.0, 0.0],
                      hydrate=[-200000.0, 200000.0],
                      add_domain=[500000.0, 1000000.0, 0],
                      ):
        with open(filename2, 'w') as f2:
            with open(self.file, 'r+') as f1:
                f_iter = iter(f1.readline, '')
                for f_single in f_iter:
                    if f_single.startswith('ELEME'):
                        f2.write('ELEME\n')
                        while True:
                            line = f_iter.__next__()
                            if line.strip():

                                x, y, z = self.location(line)
                                sequence = bisect_left(self.xposition, x)
                                kx = self.k[sequence]
                                sub_top = self.ztop[sequence]
                                c = z - kx * x
                                delta_z = sub_top - c

                                if bound[0] <= delta_z <= bound[1]:
                                    line2 = self.set_material(line, 3)
                                    f2.write(line2)

                                elif add_domain[0] <= delta_z <= add_domain[1]:
                                    line2 = self.set_material(line, add_domain[2])
                                    f2.write(line2)

                                elif hydrate[0] <= delta_z <= hydrate[1]:
                                    f2.write(line)

                                else:
                                    line2 = self.set_material(line, 2)
                                    f2.write(line2)

                            else:
                                f2.write('\n')
                                break

                    else:
                        f2.write(f_single)

                mesh = self
                mesh.file = filename2

                return mesh

    def incon_data_sim(self, x, z, sequence):
        """
        simple way to generate P T data
        :param sequence: parameter that indicate the subdomain position of point
        :param x: current x value
        :param z: current z value
        :return: return computable P and T value
        """
        delta_z = self.zmax - z
        Pp = self.ini_P + 1022 * 9.8 * delta_z
        P = '%.3g' % Pp

        delta_z2 = self.ztop[sequence] - (z - self.k[sequence] * x)
        Tt = self.ini_T + delta_z2 * 0.03
        T = '%.3g' % Tt

        return P, T

    def write_incon_data_sim(self):
        """
        calculate the ini_P and ini_T. And output that data with class <'dict'>.
        :return: return two type of data that can fit different plot function
        """

        with open(self.file) as f:
            model_incon = {}
            id_zone = 1
            f_iter = iter(f.readline, '')

            for f_single in f_iter:
                if f_single.startswith('ELEME'):
                    while True:
                        line = f_iter.__next__()

                        if line.strip():
                            elname = line[0:5]
                            MA12 = line[15:20]
                            X, Y, Z = self.location(line)
                            elem_activity = line[80:82]
                            sequence_x = bisect_left(self.xposition, X)

                            if elem_activity.strip():
                                active = False
                            else:
                                active = True

                            P1, T1 = self.incon_data_sim(X, Z, sequence_x)
                            P = float(P1)
                            T = float(T1)

                            if self.is_number(MA12):
                                MAT = int(MA12)
                            else:
                                MAT = self.material[MA12]

                            sequence = MAT - 1
                            Sh = self.Shyd * self.h[sequence]
                            Sg = self.Sgas * self.g[sequence]
                            poro3 = self.poro1 * self.poro2[sequence]

                            EOS = self.eos[MAT]

                            if MAT == 3.0:
                                P = None
                                T = None
                                Sh = None
                                Sg = None
                                po = None
                                chem = None
                                Sa = None
                            else:
                                po = poro3
                                chem = self.c
                                Sa = 1 - Sh - Sg

                            format_zone = [id_zone,
                                           elname,
                                           MAT,
                                           X,
                                           Y,
                                           Z,
                                           P,
                                           T,
                                           po,
                                           Sh,
                                           EOS,
                                           Sg,
                                           chem,
                                           Sa,
                                           active
                                           ]

                            zone = {elname: format_zone}
                            model_incon.update(zone)
                            id_zone += 1

                        else:
                            break

        return model_incon

    def get_column_row(self):
        """
        :return: total: the total number of element in the MESH
                 row_number: the number of row
                 column_number: the number of column
                 condition: condition of distribution of MESH
        """
        l_x = []
        l_z = []

        with open(self.file) as f:
            f_iter = iter(f.readline, '')

            for f_single in f_iter:
                if f_single.startswith('ELEME'):
                    while True:
                        line = f_iter.__next__()

                        if line.strip():
                            location_x = float(line[50:60])
                            location_z = float(line[70:80])
                            l_x.append(location_x)
                            l_z.append(location_z)

                        else:
                            break

        npx = np.array(l_x)
        npz = np.array(l_z)

        x1 = npx[1]
        x2 = npx[2]

        if x1 == x2:
            condition = 'column_first'

        else:
            condition = 'row_first'

        row_number = len(np.unique(npz))
        column_number = len(np.unique(npx))
        total_mesh = len(npx)

        return total_mesh, row_number, column_number, condition

    @staticmethod
    def transform_list(list_ini):
        x_c = []
        z_c = []
        p = []
        t = []
        sh = []
        sg = []
        sa = []
        chem = []
        mat = []
        for list_single in list_ini:
            mat.append(list_single[2])
            x_c.append(list_single[3])
            z_c.append(list_single[5])
            p.append(list_single[6])
            t.append(list_single[7])
            sh.append(list_single[9])
            sg.append(list_single[11])
            sa.append(list_single[13])
            chem.append(list_single[12])
        return x_c, z_c, p, t, sh, sg, chem, sa, mat

    @staticmethod
    def transform_dict(dict_ini):
        x_c = []
        z_c = []
        p = []
        t = []
        sh = []
        sg = []
        sa = []
        chem = []
        for dict_single in dict_ini:
            x_c.append(dict_ini[dict_single][3])
            z_c.append(dict_ini[dict_single][5])
            p.append(dict_ini[dict_single][6])
            t.append(dict_ini[dict_single][7])
            sh.append(dict_ini[dict_single][9])
            sg.append(dict_ini[dict_single][11])
            sa.append(dict_ini[dict_single][13])
            chem.append(dict_ini[dict_single][12])
        return x_c, z_c, p, t, sh, sg, chem, sa

    def plot_data_2d(self, dict_ini=None):
        dict_plot = dict_ini if dict_ini else self.write_incon_data_sim()

        list_rearrange = sorted(dict_plot.values(), key=(lambda d: [d[3], d[5]]))
        x, z, p, t, sh, sg, chem, sa, mat = self.transform_list(list_rearrange)

        npx1 = np.array(x)
        npz1 = np.array(z)

        x_number = len(np.unique(npx1))
        z_number = len(np.unique(npz1))

        npx = np.array(x).reshape(x_number, z_number).transpose()
        npz = np.array(z).reshape(x_number, z_number).transpose()
        p_xz = np.array(p).reshape(x_number, z_number).transpose()
        t_xz = np.array(t).reshape(x_number, z_number).transpose()
        shyd_xz = np.array(sh).reshape(x_number, z_number).transpose()
        sgas_xz = np.array(sg).reshape(x_number, z_number).transpose()
        saqu_xz = np.array(sa).reshape(x_number, z_number).transpose()
        mat_xz = np.array(mat).reshape(x_number, z_number).transpose()
        chem_xz = np.array(chem).reshape(x_number, z_number).transpose()

        dict_plot_data = {'p': p_xz,
                          't': t_xz,
                          'Shyd': shyd_xz,
                          'Sgas': sgas_xz,
                          'Saqu': saqu_xz,
                          'x_pos': npx,
                          'z_pos': npz,
                          'mat': mat_xz,
                          'chem': chem_xz,
                          }

        return dict_plot_data

    def ini_plot(self, dict_ini=None, chemical=False):

        dirname = os.path.split(self.file)[0]

        dict_plot_data = self.plot_data_2d(dict_ini)

        npx = dict_plot_data['x_pos']
        npz = dict_plot_data['z_pos']
        p_xz = dict_plot_data['p']
        t_xz = dict_plot_data['t']
        shyd_xz = dict_plot_data['Shyd']
        sgas_xz = dict_plot_data['Sgas']
        saqu_xz = dict_plot_data['Saqu']
        mat_xz = dict_plot_data['mat']
        chem_xz = dict_plot_data['chem']

        pic_name1 = 'ini_P.jpg'
        pic_path1 = os.path.join(dirname, pic_name1)
        plt.subplots(figsize=(14, 6))
        plt.contourf(npx,
                     npz,
                     p_xz,
                     40,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.title('Initial pressure')
        plt.savefig(pic_path1, dpi=1000)
        plt.show()

        pic_name2 = 'ini_T.jpg'
        pic_path2 = os.path.join(dirname, pic_name2)
        plt.subplots(figsize=(14, 6))
        plt.contourf(npx,
                     npz,
                     t_xz,
                     40,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.title('Initial temperature')
        plt.savefig(pic_path2, dpi=1000)
        plt.show()

        pic_name3 = 'ini_Shyd.jpg'
        pic_path3 = os.path.join(dirname, pic_name3)
        plt.subplots(figsize=(14, 6))
        plt.contourf(npx,
                     npz,
                     shyd_xz,
                     400,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.title('Initial hydrate saturation')
        plt.savefig(pic_path3, dpi=1000)
        plt.show()

        pic_name4 = 'ini_Sgas.jpg'
        pic_path4 = os.path.join(dirname, pic_name4)
        plt.subplots(figsize=(14, 6))
        plt.contourf(npx,
                     npz,
                     sgas_xz,
                     400,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.title('Initial gas saturation')
        plt.savefig(pic_path4, dpi=1000)
        plt.show()

        pic_name5 = 'ini_Saqu.jpg'
        pic_path5 = os.path.join(dirname, pic_name5)
        plt.subplots(figsize=(14, 6))
        plt.contourf(npx,
                     npz,
                     saqu_xz,
                     400,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.title('Initial aqua saturation')
        plt.savefig(pic_path5, dpi=1000)
        plt.show()

        pic_name6 = 'ini_mat.jpg'
        pic_path6 = os.path.join(dirname, pic_name6)
        plt.subplots(figsize=(14, 6))
        plt.contour(npx,
                    npz,
                    mat_xz,
                    3,
                    colors='k',
                    )
        plt.xticks([])
        plt.yticks([])
        plt.savefig(pic_path6, dpi=1000)
        plt.show()

        if chemical:
            pic_name7 = 'ini_chem.jpg'
            pic_path7 = os.path.join(dirname, pic_name7)
            plt.subplots(figsize=(14, 6))
            plt.contourf(npx,
                         npz,
                         chem_xz,
                         400,
                         cmap='jet',
                         )
            plt.colorbar()
            plt.xlabel('Length, x(m)')
            plt.ylabel('Elevation below top of model, z(m)')
            plt.title('Initial salinity')
            plt.savefig(pic_path7, dpi=1000)
            plt.show()

    @staticmethod
    def condition(temperature, pressure, salinity):

        pressure = pressure / 1e6
        temperature = temperature + 273.0
        intercept = 0.00059 * (salinity / 3.0) ** 2 + 0.00253 * (salinity / 3.0) + 1
        L1 = math.log(pressure / 2.23) + 35.0 * (273 / temperature) - 35 * intercept
        L2 = math.log(pressure / 2.23) + 7.5 * (273 / temperature) - 7.5 * intercept

        if L1 > 0 and L2 > 0:
            condition = 'AqH'

        elif L1 == 0 and L2 >= 0:
            condition = 'AGH'

        elif L1 >= 0 and L2 == 0:
            condition = 'AGH'

        else:
            condition = 'AqG'

        return condition

    def check_thc(self,
                  dict_ini=None):

        dict_data = dict_ini if dict_ini else self.write_incon_data_sim()[0]
        wrong_eleme = []
        dict_wrong = {}

        for data_key in dict_data:
            data_value = dict_data[data_key]
            data_p = data_value[6]
            data_t = data_value[7]
            data_c = data_value[12]
            data_sh = data_value[9]
            data_condition = data_value[10]
            data_label = data_value[1]

            if data_sh != 0:
                condition = self.condition(data_t, data_p, data_c)

                if condition != data_condition:
                    wrong_eleme.append(data_label)

        if not wrong_eleme:
            print('every element is ok')
        else:
            print('wrong element is :' + str(wrong_eleme))

        for i in wrong_eleme:
            dict_wrong.update(dict_data[i])

    def get_save(self,
                 SAVE_file='SAVE',
                 dict_ini=False,
                 chemical=True):
        """
        get data from SAVE file, it needs the first type of dict from self.write_incon_data_sim or
         self.write_noincon_data. if data is from write_incon_data_sim, we don't need to input dict because we have
         invoked it in this function. if data is not from ..., we need to get the data and input it
        :param SAVE_file: SAVE file
        :param dict_ini: dict that we need in this function
        :param chemical: if we need to consider the influence of chemical
        :return: return a new dict with data from SAVE file
        """
        if not dict_ini:
            dict_ini = self.write_incon_data_sim()[0]

        dirname = os.path.split(self.file)[0]

        save_f = os.path.join(dirname, SAVE_file) if SAVE_file == 'SAVE' else SAVE_file

        dict_new = dict_ini

        with open(save_f) as f:
            f_iter = iter(f.readline, '')
            for f_single in f_iter:
                if f_single.startswith('INCON'):
                    while True:
                        line = f_iter.__next__()

                        if line.startswith(':::'):
                            break

                        else:
                            label = line[0:5]
                            state = line[32:35]
                            line2 = f_iter.__next__()

                            if chemical:
                                if state == 'Aqu':
                                    p = '%.3g' % float(line2[:20])
                                    sg = '%.3g' % 0
                                    sh = '%.3g' % float(line2[40:60])
                                    c = '%.3g' % float(line2[60:80])
                                    t = '%.3g' % float(line2[80:100])

                                elif state == 'AqH':
                                    p = '%.3g' % float(line2[:20])
                                    sg = '%.3g' % 0
                                    sh = '%.3g' % (1 - float(line2[20:40]))
                                    c = '%.3g' % float(line2[60:80])
                                    t = '%.3g' % float(line2[80:100])

                                elif state == 'AqG':
                                    p = '%.3g' % float(line2[:20])
                                    sg = '%.3g' % (1 - float(line2[20:40]) - float(line2[40:60]))
                                    sh = '%.3g' % float(line2[40:60])
                                    c = '%.3g' % float(line2[60:80])
                                    t = '%.3g' % float(line2[80:100])

                                elif state == 'AGH':
                                    p = '%.3g' % float(line2[:20])
                                    sg = '%.3g' % float(line2[40:60])
                                    sh = '%.3g' % (1 - float(line2[40:60]) - float(line2[20:40]))
                                    c = '%.3g' % float(line2[60:80])
                                    t = '%.3g' % float(line2[80:100])

                            else:
                                if state == 'Aqu':
                                    p = '%.3g' % float(line2[:20])
                                    sg = '%.3g' % 0
                                    sh = '%.3g' % float(line2[40:60])
                                    t = '%.3g' % float(line2[60:80])

                                elif state == 'AqH':
                                    p = '%.3g' % float(line2[:20])
                                    sg = '%.3g' % 0
                                    sh = '%.3g' % (1 - float(line2[20:40]))
                                    t = '%.3g' % float(line2[60:80])

                                elif state == 'AqG':
                                    p = '%.3g' % float(line2[:20])
                                    sg = '%.3g' % (1 - float(line2[20:40]) - float(line2[40:60]))
                                    sh = '%.3g' % float(line2[40:60])
                                    t = '%.3g' % float(line2[60:80])

                                elif state == 'AGH':
                                    p = '%.3g' % float(line2[:20])
                                    sg = '%.3g' % float(line2[40:60])
                                    sh = '%.3g' % (1 - float(line2[40:60]) - float(line2[20:40]))
                                    t = '%.3g' % float(line2[60:80])

                            if dict_new[label][14]:

                                dict_new[label][6] = float(p)
                                dict_new[label][7] = float(t)
                                dict_new[label][9] = float(sh)
                                dict_new[label][11] = float(sg)
                                dict_new[label][13] = float(1 - float(sh) - float(sg))
                                if chemical:
                                    dict_new[label][12] = float(c)
                            else:
                                dict_new[label][6] = None
                                dict_new[label][7] = None
                                dict_new[label][9] = None
                                dict_new[label][11] = None
                                dict_new[label][13] = None
                                if chemical:
                                    dict_new[label][12] = None

        return dict_new

    def save_plot(self,
                  SAVE_file,
                  dict_save=False,
                  chemical=True):
        """
        plot data of result. data is got from self.save_plot, if data is False, invoke this function.

        :param SAVE_file: SAVE file
        :param dict_save: class <'dict'> that contain the data of save file. it is got from function self.get_save
        :param chemical: if we consider the influence of chemical
        :return: no return
        """
        if not dict_save:
            dict_save = self.get_save(SAVE_file, False, chemical=chemical)

        dirname = os.path.split(self.file)[0]

        dict_plot_data = self.plot_data_2d(dict_save)
        npx = dict_plot_data['x_pos']
        npz = dict_plot_data['z_pos']
        p_xz = dict_plot_data['p']
        t_xz = dict_plot_data['t']
        shyd_xz = dict_plot_data['Shyd']
        sgas_xz = dict_plot_data['Sgas']
        saqu_xz = dict_plot_data['Saqu']
        chem_xz = dict_plot_data['chem']

        pic_name1 = 'Save_P.jpg'
        pic_path1 = os.path.join(dirname, pic_name1)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     p_xz,
                     40,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.title('Saved pressure')
        plt.savefig(pic_path1, dpi=1000)
        plt.show()

        pic_name2 = 'Save_T.jpg'
        pic_path2 = os.path.join(dirname, pic_name2)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     t_xz,
                     40,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.title('Saved temperature')
        plt.savefig(pic_path2, dpi=1000)
        plt.show()

        pic_name3 = 'Save_Shyd.jpg'
        pic_path3 = os.path.join(dirname, pic_name3)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     shyd_xz,
                     20,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.title('Saved hydrate saturation')
        plt.savefig(pic_path3, dpi=1000)
        plt.show()

        pic_name4 = 'Save_Sgas.jpg'
        pic_path4 = os.path.join(dirname, pic_name4)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     sgas_xz,
                     400,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.title('Saved gas saturation')
        plt.savefig(pic_path4, dpi=1000)
        plt.show()

        pic_name5 = 'Save_Saqu.jpg'
        pic_path5 = os.path.join(dirname, pic_name5)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     saqu_xz,
                     20,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.title('Saved aqua saturation')
        plt.savefig(pic_path5, dpi=1000)
        plt.show()

        if chemical:
            pic_name6 = 'Save_chem.jpg'
            pic_path6 = os.path.join(dirname, pic_name6)
            plt.subplots(figsize=(9, 6))
            plt.contourf(npx,
                         npz,
                         chem_xz,
                         20,
                         cmap='jet',
                         )
            plt.colorbar()
            plt.xlabel('Length, x(m)')
            plt.ylabel('Elevation below top of model, z(m)')
            plt.title('Saved salinity')
            plt.savefig(pic_path6, dpi=1000)
            plt.show()

    def elem_find(self, x_value, z_value, y_value=None, dict_ini=None):

        label = None

        dict_infor = dict_ini if dict_ini else self.write_incon_data_sim()

        for k, v in dict_infor.items():

            if v[3] == x_value:

                if v[5] == z_value:

                    if y_value:

                        if v[4] == y_value:
                            label = k

                    else:
                        label = k

        try:
            label

        except NameError:
            raise ValueError("There's no elem with indicated xyz_value")

        return label

    def local_plot(self, z_range, dict_ini=None, x_range=None, y_range=None):

        dict_global = dict_ini if dict_ini else self.write_incon_data_sim()[0]
        dict_local = {}

        dirname = os.path.split(self.file)[0]

        for k, v in dict_global.items():
            if z_range[0] <= v[5] <= z_range[1]:
                if x_range:
                    if y_range:
                        if x_range[0] <= v[3] <= x_range[1] and y_range[0] <= v[4] <= y_range[1]:
                            local_key = k
                            local_value = v
                            single_local = {local_key: local_value}
                            dict_local.update(single_local)

                    else:
                        if x_range[0] <= v[3] <= x_range[1]:
                            local_key = k
                            local_value = v
                            single_local = {local_key: local_value}
                            dict_local.update(single_local)
                else:

                    local_key = k
                    local_value = v
                    single_local = {local_key: local_value}
                    dict_local.update(single_local)

        local_plot_data = self.plot_data_2d(dict_local)
        npx = local_plot_data['x_pos']
        npz = local_plot_data['z_pos']
        p_xz = local_plot_data['p']
        t_xz = local_plot_data['t']
        shyd_xz = local_plot_data['Shyd']
        sgas_xz = local_plot_data['Sgas']
        saqu_xz = local_plot_data['Saqu']
        chem_xz = local_plot_data['chem']

        pic_name1 = 'local_P.jpg'
        pic_path1 = os.path.join(dirname, pic_name1)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     p_xz,
                     40,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.savefig(pic_path1, dpi=1000)
        plt.show()

        pic_name2 = 'local_T.jpg'
        pic_path2 = os.path.join(dirname, pic_name2)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     t_xz,
                     40,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.savefig(pic_path2, dpi=1000)
        plt.show()

        pic_name3 = 'local_Shyd.jpg'
        pic_path3 = os.path.join(dirname, pic_name3)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     shyd_xz,
                     400,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.savefig(pic_path3, dpi=1000)
        plt.show()

        pic_name4 = 'local_Sgas.jpg'
        pic_path4 = os.path.join(dirname, pic_name4)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     sgas_xz,
                     400,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.savefig(pic_path4, dpi=1000)
        plt.show()

        pic_name5 = 'local_Saqu.jpg'
        pic_path5 = os.path.join(dirname, pic_name5)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     saqu_xz,
                     400,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.savefig(pic_path5, dpi=1000)
        plt.show()

        pic_name6 = 'local_chem.jpg'
        pic_path6 = os.path.join(dirname, pic_name6)
        plt.subplots(figsize=(9, 6))
        plt.contourf(npx,
                     npz,
                     chem_xz,
                     400,
                     cmap='jet',
                     )
        plt.colorbar()
        plt.xlabel('Length, x(m)')
        plt.ylabel('Elevation below top of model, z(m)')
        plt.savefig(pic_path6, dpi=1000)
        plt.show()

    def read_time_series(self, merge=True):
        dirname = os.path.split(self.file)[0]

        for path in os.listdir(dirname):
            if r'Time_Series' in path:
                path_ts = path
                abs_path = os.path.join(dirname, path_ts)

                with open(abs_path) as f:
                    first_line = f.readline()

                well_data = pd.read_fwf(abs_path, skiprows=1)
                well_data2 = well_data[['Time:days', 'M_A:kg', 'M_CH4G:kg', 'Q_A:kg/s', 'Q_CH4G:kg/s']]
                time = np.array(well_data2['Time:days'])
                M_A = np.array(well_data2['M_A:kg'])
                Q_A = np.array(well_data2['Q_A:kg/s']) * 86400 / 1022
                M_G = np.array(well_data2['M_CH4G:kg'])
                Q_G = np.array(well_data2['Q_CH4G:kg/s']) * 86400 / 0.716
                V_A = M_A / 1022
                V_G = M_G / 0.716

                if merge:

                    fig, ax = plt.subplots()
                    fig.set_size_inches(8, 4.8)
                    pic_name1 = 'Qa&Va.jpg'
                    pic_path1 = os.path.join(dirname, pic_name1)
                    ax.plot(time, V_A, color='blue', label='V_A', linewidth=2.0)
                    ax.plot(time, Q_A, color='red', label='Q_A', linewidth=2.0)
                    ax.set_title('water')
                    ax.grid(which='both', linestyle='--', linewidth=1.5)
                    ax.set_xlabel('Time:days', fontsize=15)
                    ax.set_ylabel('V_A(m³) & Q_A(m³/day)', fontsize=15)
                    ax.spines['bottom'].set_linewidth('2.0')
                    ax.spines['top'].set_linewidth('2.0')
                    ax.spines['right'].set_linewidth('2.0')
                    ax.spines['left'].set_linewidth('2.0')
                    ax.legend()
                    plt.savefig(pic_path1, dpi=1000)
                    plt.show()

                    fig, ax = plt.subplots()
                    fig.set_size_inches(8, 4.8)
                    pic_name2 = 'Qg&Vg.jpg'
                    pic_path2 = os.path.join(dirname, pic_name2)
                    ax.plot(time, V_G, color='blue', label='V_G', linewidth=2.0)
                    ax.plot(time, Q_G, color='red', label='Q_G', linewidth=2.0)
                    ax.set_title('gas')
                    ax.grid(which='major', linestyle='--', linewidth=1.0)
                    ax.set_xlabel('Time:days', fontsize=15)
                    ax.set_ylabel('V_G(m³) & Q_G(m³/day)', fontsize=15)
                    ax.spines['bottom'].set_linewidth('2.0')
                    ax.spines['top'].set_linewidth('2.0')
                    ax.spines['right'].set_linewidth('2.0')
                    ax.spines['left'].set_linewidth('2.0')
                    plt.savefig(pic_path2, dpi=1000)
                    ax.legend()
                    plt.show()

                else:

                    fig, ax = plt.subplots()
                    fig.set_size_inches(7, 4.8)
                    pic_name1 = 'Qa.jpg'
                    pic_path1 = os.path.join(dirname, pic_name1)
                    ax.plot(time, Q_A, color='red', label='Q_A')
                    ax.set_xlabel('Time:days')
                    ax.set_ylabel('Q_A(m3/day)')
                    ax.legend()
                    plt.savefig(pic_path1, dpi=1000)
                    plt.show()

                    fig, ax = plt.subplots()
                    fig.set_size_inches(7, 4.8)
                    pic_name2 = 'Va.jpg'
                    pic_path2 = os.path.join(dirname, pic_name2)
                    ax.plot(time, V_A, color='red', label='V_A')
                    ax.set_xlabel('Time:days')
                    ax.set_ylabel('V_A(m3/day)')
                    ax.legend()
                    plt.savefig(pic_path2, dpi=1000)
                    plt.show()

                    fig, ax = plt.subplots()
                    fig.set_size_inches(7, 4.8)
                    pic_name3 = 'Qg.jpg'
                    pic_path3 = os.path.join(dirname, pic_name3)
                    ax.plot(time, Q_A, color='red', label='Q_G')
                    ax.set_xlabel('Time:days')
                    ax.set_ylabel('Q_G(m3/day)')
                    ax.legend()
                    plt.savefig(pic_path3, dpi=1000)
                    plt.show()

                    fig, ax = plt.subplots()
                    fig.set_size_inches(7, 4.8)
                    pic_name4 = 'Vg.jpg'
                    pic_path4 = os.path.join(dirname, pic_name4)
                    ax.plot(time, V_A, color='red', label='V_G')
                    ax.set_xlabel('Time:days')
                    ax.set_ylabel('V_G(m3/day)')
                    ax.legend()
                    plt.savefig(pic_path4, dpi=1000)
                    plt.show()
