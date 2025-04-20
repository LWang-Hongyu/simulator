#!/bin/bash

# 确保脚本在遇到错误时停止执行
set -e

#读取变量
suffix=$(head -n 1 /home/why/DNNFlow/simulation/mix/flow.txt)

#检查是否成功读取到suffix
if [ -z "$suffix" ]; then
        echo "Error: Failed to read the first line from /home/why/DNNFlow/simulation/mix/flow.txt"
        exit 1
fi

# 重新配置
echo "Restart configure..."
sudo ./waf configure

# 执行waf命令
echo "Running waf command..."
sudo ./waf --run 'scratch/third mix/config.txt' > log826.txt

scp log826.txt zdw@115.157.197.176:/home/zdw/

# 切换到analysis目录
echo "Changing to analysis directory..."
cd ../analysis/

# 执行trace_reader命令，并将输出重定向到文件
echo "Running trace_reader and redirecting output to flow_${suffix}.txt..."
sudo ./trace_reader ../simulation/mix/mix.tr > flow_${suffix}.txt

# 切换到visualization/src目录
echo "Changing to visualization/src directory..."
cd ../visualization/src/

# 执行Python脚本
echo "Running flow_visual.py..."
sudo python3 flow_visual.py -i ../../analysis/flow_${suffix}.txt

echo "Script execution completed successfully."

#把生成的文件拷贝到.2上
echo "copy png to zdw@115.157.197.176..."
cd output/flow_${suffix}/

scp flow-25000.png zdw@115.157.197.176:/home/zdw/
