import os


def get_dir_files(dir_name=None):
    """to gain the input file list

    :parameter
        dir_name:string
            directory name or file name

    :returns
        file_list: list
            file list contain all the absolute file path under the given dir
    """
    import os

    file_list = []
    if not os.path.exists(dir_name):
        return 0
    if os.path.isfile(dir_name):
        return file_list.append(dir_name)
    else:
        for strRoot, lsDir, lsFiles in os.walk(dir_name):
            for fn in lsFiles:
                file_list.append(os.path.join(strRoot, fn))
            if len(lsDir) > 0:
                for dn in lsDir:
                    file_list.append(os.path.join(strRoot, dn))
    return file_list


def is_handled(dir_name):
    for root, ls_dir, ls_files in os.walk(dir_name):
        for fi in ls_files:
            if fi.split('.')[-1] == "mat":
                return True
        return False


def get_dir(dir_name):
    n = 0
    dir_list = []
    dir_unhandle = []
    for str_root, ls_dir, ls_files in os.walk(dir_name):
        for dn in ls_dir:
            dir_list.append(os.path.join(str_root,dn))

    #a = [len(get_dir_files(item)) for item in  dir_list]
    for item in dir_list:
        files = get_dir_files(item)
        a = [len(files)]
    if len(files) < 100:
        print item, '\t', dir_list.index(item), '\t', len(files)
        n += 1
    print sorted(a)

    # print n


if __name__ == '__main__':
    get_dir("E:\\graph_matching\\fgm\\src\\asg\\fgm\\data_process_Nattr")
