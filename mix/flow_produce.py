import sys
import random

src = 0
dst = 64


def write_random_numbers_to_file(line_count, filename):
    global src, dst
    with open(filename, 'w') as file:
        file.write(str(line_count) + '\n')
        for _ in range(line_count):
            # 生成一个0到1之间的随机浮点数,表示通信时间占用的比例
            ratio = random.uniform(0.4,0.6)
            bits_per_iteration = random.randint(1e7, 1e8)
            iteration = random.randint(1e1, 1e1)
            total_bytes = iteration * bits_per_iteration
            gap_per_iteration = bits_per_iteration * 8 / 1e11 / ratio * (1 - ratio) * 1e9
            # 将随机数写入文件
            file.write(str(src) + ' ' + str(dst) + ' 3 100 ' + str(total_bytes) + ' 0 ')
            src += 1
            dst += 1
            file.write(str(bits_per_iteration) + ' ' + str(round(gap_per_iteration)) + '\n')
            #file.write(str(iteration) + ' ' + str(ratio) + '\n')


def main():
    # 检查命令行参数的数量
    if len(sys.argv) != 2:
        print("Usage: python flow_produce.py <number_of_lines>")
        sys.exit(1)

    try:
        # 将命令行参数转换为整数
        line_count = int(sys.argv[1])
        if line_count <= 0:
            raise ValueError("Number of lines must be a positive integer.")
    except ValueError as e:
        print("Error:", e)
        print("Usage: python script.py <number_of_lines>")
        sys.exit(1)

    # 指定输出文件名
    filename = "flow.txt"

    # 生成随机数并写入文件
    write_random_numbers_to_file(line_count, filename)
    print(f"Wrote {line_count} random numbers to {filename}")


if __name__ == "__main__":
    main()
