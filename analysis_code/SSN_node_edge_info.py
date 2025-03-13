# 绘制DNB的landscape
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D  


# 提取DNB基因
def read_csv(stage_id, sample_num):
    file = r"/Pancreatitis/Huaxi_Pancreatitis/lDNB_output/Max_dnb_score_module in {} for sample {}.csv".format(stage_id,
                                                                                                     sample_num)
    df = pd.read_csv(file, header=0, index_col=0)
    return df.index[0:100]


# 求至少三个样本中表达的基因
def sample_union(stage_id, sample_num):
    variable_name = [None] * sample_num
    tmp = []
    for num in range(1, sample_num + 1):
        variable_name[num - 1] = read_csv(stage_id, num)
        tmp.append(variable_name[num - 1])
        if num == 1:
            r = list(variable_name[0])
        else:
            r = list(set(variable_name[0]).union(variable_name[num - 1]))
    return r, tmp


Health_dnb, Health_dnb_lst = sample_union('Health', 5)
MAP_dnb, MAP_dnb_lst = sample_union('MAP', 56)
SAP_dnb, SAP_dnb_lst = sample_union('SAP', 12)
print('Health dnb len: ', len(Health_dnb))
print('MAP dnb len: ', len(MAP_dnb))
print('SAP dnb len: ', len(SAP_dnb))
print('Health dnb len: ', len(Health_dnb_lst))
print('MAP dnb len: ', len(MAP_dnb_lst))
print('SAP dnb len: ', len(SAP_dnb_lst))

# 绘制DNB单样本网络图
# 边界点取并集
def ssn_edge_union(stage_id, sample_num):
    variable_name = [None] * sample_num
    ssn_edge_temp = [None] * sample_num
    for num in range(1, sample_num + 1):
        variable_name[num - 1] = pd.read_csv(
            '/Pancreatitis/Huaxi_Pancreatitis/lDNB_output/SSN for ' + stage_id + ' in sample ' + str(num) + '.csv')
        ssn_edge_temp[num - 1] = variable_name[num - 1].iloc[:, 0] + '_' + variable_name[num - 1].iloc[:, 1]
        if num == 1:
            r = list(ssn_edge_temp[0])
        else:
            r = list(set(ssn_edge_temp[0]).union(ssn_edge_temp[num - 1]))
    return r


def get_SSN_edge(stage_id, sample_num):
    """
        计算目标时间点所有样本的SSN：具体操作为保留至少出现在两个样本的边
    """
    df = pd.DataFrame(columns=['node1', 'node2'])
    ssn_edge = ssn_edge_union(stage_id, sample_num)
    for edge in ssn_edge:
        df = pd.concat([df, pd.DataFrame({'node1': [edge.split('_')[0]], 'node2': [edge.split('_')[1]]})], ignore_index=True)
    df.to_csv(r"/Pancreatitis/Huaxi_Pancreatitis/lDNB_output//SSN for AP data ' + str(stage_id) + '.csv", header=True, index=False)
    print('ssn edge num: ', len(ssn_edge))
    return ssn_edge


sample_num_lst = [5, 56, 12]  # 每个stage的样本个数
stage_lst = ['Health', 'MAP', 'SAP']
total_ssn = pd.DataFrame(columns=['node1', 'n1_isTP1_dnb', 'n1_isTP2_dnb', 'node2', 'n2_isTP1_dnb', 'n2_isTP2_dnb',
                                  'Health', 'MAP', 'SAP'])
for i in range(3):
    print('当前stage为：', stage_lst[i])
    temp_ssn_edge = get_SSN_edge(stage_lst[i], sample_num_lst[i])
    for edge in temp_ssn_edge:
        if edge not in total_ssn.index:
            total_ssn.loc[edge, ['node1', 'node2', 'stage' + str(i + 1)]] = [edge.split('_')[0], edge.split('_')[1], 1]
        else:
            total_ssn.loc[edge, ['stage' + str(i + 1)]] = 1

print('ssn background edge num: ', total_ssn.shape[0])

total_ssn.to_csv(r"/Pancreatitis/Huaxi_Pancreatitis/lDNB_output/total_ssn_for_AP data.csv")
ssn_edge_info = total_ssn.drop(columns=['n1_isTP1_dnb', 'n1_isTP2_dnb', 'n2_isTP1_dnb', 'n2_isTP2_dnb'])
ssn_edge_info.fillna(0).to_csv(r"/Pancreatitis/Huaxi_Pancreatitis/lDNB_output/ssn_edge_info.csv", index=False)

# 设置background network node的属性
ssn_node_info = pd.DataFrame(columns=['node_id', 'node_type'])
ssn_node_info.loc[:, 'node_id'] = list(set(total_ssn.loc[:, 'node1']) | set(total_ssn.loc[:, 'node2']))
ssn_node_info.loc[:, ['node_type']] = 0

def stage_score_union(stage_id, sample_num):
    variable_name = [None] * sample_num
    for num in range(1, sample_num + 1):
        file = r"/Pancreatitis/Huaxi_Pancreatitis/lDNB_output/Max_dnb_score_module in {} for sample {}.csv".format(stage_id, num)
        variable_name[num - 1] = pd.read_csv(file, header=0, index_col=0)
        if num == 1:
            r = list(variable_name[0].index)
        else:
            r = list(set(variable_name[0].index).union(variable_name[num - 1].index))
    return r

def stage_dnb_score(stage_id, sample_num):
    """
        计算目标时间点DNB基因的score值，这里仅计算至少在两个样本中出现的DNB基因
    """
    idx = stage_score_union(stage_id, sample_num)
    df = pd.DataFrame(index=idx, columns=[x for x in range(1, sample_num + 1)])
    for num in range(1, sample_num + 1):
        file = r"/Pancreatitis/Huaxi_Pancreatitis/lDNB_output/Max_dnb_score_module in {} for sample {}.csv".format(stage_id, num)
        variable_temp = pd.read_csv(file, header=0, index_col=0)
        idx1 = variable_temp.index
        df.loc[idx, num] = variable_temp.loc[idx1, 'score']

    df.loc[:, 'mean_score'] = df.mean(axis=1)

    return df


# 计算每个时间解阶段的dnb分值
sample_num_lst = [5, 56, 12]  # 每个stage的样本个数
for i in range(3):
    print('当前stage为：', stage_lst[i])
    temp_stage_score = stage_dnb_score(stage_lst[i], sample_num_lst[i])
    in_set = set(ssn_node_info['node_id']) & set(temp_stage_score.index)
    ssn_node_info.loc[ssn_node_info.loc[:, 'node_id'].isin(temp_stage_score.index), 'T' + str(i + 1)] = \
        temp_stage_score.loc[in_set, 'mean_score'].values

ssn_node_info = ssn_node_info.fillna(0)
ssn_node_info.to_csv(r"/Pancreatitis/Huaxi_Pancreatitis/lDNB_output/ssn_node_info.csv", index=False)

