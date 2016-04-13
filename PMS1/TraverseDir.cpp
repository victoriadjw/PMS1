#include"TraverseDir.h"


CBrowseDir::CBrowseDir(char *of, char *f)
{
	//�õ�ǰĿ¼��ʼ��m_szInitDir
	ofs.open(output_file, ios::app | ios::out);

	_getcwd(m_szInitDir, _MAX_PATH);

	//���Ŀ¼�����һ����ĸ����'\',����������һ��'\'
	int len = strlen(m_szInitDir);
	if (m_szInitDir[len - 1] != '\\')
		strcat(m_szInitDir, "\\");
}
CBrowseDir::CBrowseDir()
{
	//�õ�ǰĿ¼��ʼ��m_szInitDir

	_getcwd(m_szInitDir, _MAX_PATH);

	//���Ŀ¼�����һ����ĸ����'\',����������һ��'\'
	int len = strlen(m_szInitDir);
	if (m_szInitDir[len - 1] != '\\')
		strcat(m_szInitDir, "\\");
}

bool CBrowseDir::SetInitDir(const char *dir)
{
	//�Ȱ�dirת��Ϊ����·��
	if (_fullpath(m_szInitDir, dir, _MAX_PATH) == NULL)
		return false;

	//�ж�Ŀ¼�Ƿ����
	if (_chdir(m_szInitDir) != 0)
		return false;

	//���Ŀ¼�����һ����ĸ����'\',����������һ��'\'
	int len = strlen(m_szInitDir);
	if (m_szInitDir[len - 1] != '\\')
		strcat(m_szInitDir, "\\");

	return true;
}

bool CBrowseDir::BeginBrowse(const char *filespec)
{
	ProcessDir(m_szInitDir, NULL);
	return BrowseDir(m_szInitDir, filespec);
}

bool CBrowseDir::BrowseDir(const char *dir, const char *filespec)
{
	_chdir(dir);

	//���Ȳ���dir�з���Ҫ����ļ�
	long hFile;
	_finddata_t fileinfo;
	if ((hFile = _findfirst(filespec, &fileinfo)) != -1)
	{
		do
		{
			//����ǲ���Ŀ¼
			//�������,����д���
			if (!(fileinfo.attrib & _A_SUBDIR))
			{
				char filename[_MAX_PATH];
				//strcpy(filename, dir);
				strcpy(filename, "");
				strcat(filename, fileinfo.name);
				cout << filename << endl;
				//ReadFile(fileinfo.name);	
				//vec_all_file.push_back(filename);
				vec_all_file.push_back(filename);
				if (!ProcessFile(filename))
					return false;
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
	//����dir�е���Ŀ¼
	//��Ϊ�ڴ���dir�е��ļ�ʱ���������ProcessFile�п��ܸı���
	//��ǰĿ¼����˻�Ҫ�������õ�ǰĿ¼Ϊdir��
	//ִ�й�_findfirst�󣬿���ϵͳ��¼���������Ϣ����˸ı�Ŀ¼
	//��_findnextû��Ӱ�졣
	_chdir(dir);
	if ((hFile = _findfirst("*.*", &fileinfo)) != -1)
	{
		do
		{
			//����ǲ���Ŀ¼
			//�����,�ټ���ǲ��� . �� .. 
			//�������,���е���
			if ((fileinfo.attrib & _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp
					(fileinfo.name, "..") != 0)
				{
					char subdir[_MAX_PATH];
					strcpy(subdir, dir);
					strcat(subdir, fileinfo.name);
					strcat(subdir, "\\");
					ProcessDir(subdir, dir);
					if (!BrowseDir(subdir, filespec))
						return false;
				}
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
	return true;
}
void CBrowseDir::ReadFile(const char *fn)
{
	string strfn = fn;
	if (strfn.find(filter) == string::npos)
		return;
	ifstream ifs(fn, ios::_Nocreate);
	if (!ifs)
	{
		cerr << "file: " << fn << " open error";
		return;
	}
	string line;
	int vs, run_times, num_iter, num_search, imp[3];
	double t, total_t;
	bool vs_once = true;// for the bug that the second VS is 0
	while (getline(ifs, line))
	{
		cout << line << endl;
		char *ch;
		const int len = line.length();
		ch = new char[len + 1];
		strcpy(ch, line.c_str());

		const char *split = " \t\f\r\v\n";
		if (line.find("VS=") != string::npos/*&&vs_once*/)
		{
			int split_cnt = 0;
			for (char *p = strtok(ch, split); p != NULL; p = strtok(NULL, split), split_cnt++)
			{
				//cout << p << " " << split_cnt << endl;
				if (split_cnt == 1)
					vs = atoi(p);
				if (split_cnt == 4)
					t = atof(p);
				if (split_cnt == 7)
					run_times = atoi(p);
			}
			vs_once = false;
		}
		if (line.find("ITERATION") != string::npos)
		{
			int split_cnt = 0;
			for (char *p = strtok(ch, split); p != NULL; p = strtok(NULL, split), split_cnt++)
			{
				if (split_cnt == 1)
					num_iter = atoi(p);
			}
		}
		if (line.find("SEARCH") != string::npos)
		{
			int split_cnt = 0;
			for (char *p = strtok(ch, split); p != NULL; p = strtok(NULL, split), split_cnt++)
			{
				if (split_cnt == 1)
					num_search = atoi(p);
			}
		}
		if (line.find("USE") != string::npos)
		{
			int split_cnt = 0;
			for (char *p = strtok(ch, split); p != NULL; p = strtok(NULL, split), split_cnt++)
			{
				if (split_cnt == 1)
					total_t = atof(p);
			}
		}
		if (line.find("improve1") != string::npos)
		{
			int split_cnt = 0;
			for (char *p = strtok(ch, split); p != NULL; p = strtok(NULL, split), split_cnt++)
			{
				if (split_cnt == 2)
					imp[0] = atof(p);
				if (split_cnt == 5)
					imp[1] = atof(p);
				if (split_cnt == 9)
					imp[2] = atof(p);
			}
			break;
		}
	}

	cout << fn << "\t"
		<< vs << "\t"
		<< t << "\t"
		<< run_times << "\t"
		<< num_iter << "\t"
		<< num_search << "\t"
		<< total_t << "\t"
		<< imp[0] << "\t"
		<< imp[1] << "\t"
		<< imp[2] << endl;

	ofs << fn << "\t"
		<< vs << "\t"
		<< t << "\t"
		<< run_times << "\t"
		<< num_iter << "\t"
		<< num_search << "\t"
		<< total_t << "\t"
		<< imp[0] << "\t"
		<< imp[1] << "\t"
		<< imp[2] << endl;
	ifs.close();
}
bool CBrowseDir::ProcessFile(const char *filename)
{
	return true;
}

void CBrowseDir::ProcessDir(const char
	*currentdir, const char *parentdir)
{
}

