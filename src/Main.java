import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

// 基因类
class Gene {
    public double fitness; // 适应度
    public Integer[] chromosome; // 染色体

    // 构造函数
    public Gene() {
        fitness = 0;
    }

    // 带参数的构造函数
    public Gene(double fitness) {
        this.fitness = fitness;
    }
}

// 结果类
class Result {
    public int fulfillTime = 0; // 完成时间
    public int[] machineWorkTime = new int[1024]; // 机器工作时间
    public int[] processIds = new int[1024]; // 过程ID
    public int[][] endTime = new int[1024][1024]; // 结束时间
    public int[][] startTime = new int[1024][1024]; // 开始时间
}

// 遗传算法类
class GeneticAlgorithm {
    private final int population_num = 60; // 种群数量
    private final double rate = 0.05; // 变异概率
    private int job_num; // 工作数量
    private int machine_num; // 机器数量
    private int process_num; // 过程数量
    private int chromosome_Size; // 染色体大小
    private int[][] machine_Matrix = new int[1024][1024]; // 机器矩阵
    private int[][] time_Matrix = new int[1024][1024]; // 时间矩阵
    private int[][] process_Matrix = new int[1024][1024]; // 过程矩阵
    private Set<Gene> geneSet = new HashSet<>(); // 基因集合
    private Random random = new Random(); // 随机数生成器

    // 构造函数
    public GeneticAlgorithm(int jobNumber, int machineNumber) {
        this.job_num = jobNumber;
        this.machine_num = machineNumber;
        this.machine_Matrix = new int[jobNumber][machineNumber];
        this.process_Matrix = new int[jobNumber][machineNumber];
        for (int i = 0; i < jobNumber; i++) {
            Arrays.fill(this.machine_Matrix[i], -1); // 将每行的元素初始化为-1
            Arrays.fill(this.process_Matrix[i], -1);
        }
    }


    // 创建从0到n-1的整数列表
    private List<Integer> makeList(int n) {
        return IntStream.range(0, n)
                .boxed()
                .collect(Collectors.toList());
    }

    //流式操作
// 从数组中过滤掉特定值
    private Integer[] filterArray(Integer[] arr, int filterVal) {
        return Arrays.stream(arr)
                .filter(val -> val != filterVal)
                .toArray(Integer[]::new);
    }

    // 从数组中获取子数组
    // 从数组中获取子数组
    private List<Integer> subArray(Integer[] arr, int start, int end) {
        int actualStart = Math.max(start, 0); // 确保开始索引不小于0
        int actualEnd = Math.min(end, arr.length); // 确保结束索引不超出数组范围
        if (actualStart >= actualEnd) {
            return new ArrayList<>(); // 如果开始索引大于等于结束索引，则返回空列表
        }
        return Arrays.stream(arr, actualStart, actualEnd)
                .collect(Collectors.toList());
    }


    // 初始化种群
    public void initialPopulation() {
        // 使用流式操作生成种群
        geneSet = IntStream.range(0, population_num) // 生成一个范围为0到population_num-1的整数流
                .mapToObj(i -> {
                    Gene g = new Gene(); // 创建一个新的基因对象
                    int size = job_num * machine_num; // 计算染色体数组的大小
                    List<Integer> indexList = makeList(size); // 生成一个包含索引的列表
                    Integer[] chromosome = new Integer[size]; // 创建一个大小为size的整数数组
                    Arrays.fill(chromosome, -1); // 将染色体数组初始化为-1

                    // 填充染色体数组
                    IntStream.range(0, job_num)
                            .forEach(j -> IntStream.range(0, machine_num)
                                    .forEach(k -> {
                                        int index = random.nextInt(indexList.size()); // 从索引列表中随机选择一个索引
                                        int val = indexList.remove(index); // 从索引列表中移除选择的索引
                                        if (process_Matrix[j][k] != -1) {
                                            chromosome[val] = j; // 根据作业分配值填充染色体
                                        }
                                    }));

                    g.chromosome = filterArray(chromosome, -1); // 过滤掉值为-1的元素
                    g.fitness = calculateFitness(g).fulfillTime; // 计算基因的适应度
                    return g; // 返回生成的基因对象
                })
                .collect(Collectors.toSet()); // 将生成的基因对象收集到一个集合中
    }


    // 计算适应度
    public Result calculateFitness(Gene g) {
        Result result = new Result(); // 创建一个用于存储结果的对象

        // 遍历基因的染色体数组
        for (int jobId : g.chromosome) {
            int processId = result.processIds[jobId]; // 获取作业的处理过程ID
            int machineId = machine_Matrix[jobId][processId]; // 获取作业在该处理过程上的机器ID
            int time = time_Matrix[jobId][processId]; // 获取作业在该处理过程上的加工时间

            updateStartTimes(result, jobId, processId, machineId); // 更新作业的开始时间
            updateMachineWorkTime(result, jobId, processId, machineId, time); // 更新机器的工作时间
            updateEndTime(result, jobId, processId, machineId); // 更新作业的结束时间

            // 更新最终完成时间为当前机器的工作时间和历史最大完成时间的较大值
            result.fulfillTime = Math.max(result.fulfillTime, result.machineWorkTime[machineId]);
        }

        return result; // 返回计算得到的结果
    }

    // 更新作业的开始时间
    private void updateStartTimes(Result result, int jobId, int processId, int machineId) {
        if (processId == 0) {
            result.startTime[jobId][processId] = result.machineWorkTime[machineId];
        } else {
            result.startTime[jobId][processId] = Math.max(result.endTime[jobId][processId - 1], result.machineWorkTime[machineId]);
        }
    }

    // 更新机器的工作时间
    private void updateMachineWorkTime(Result result, int jobId, int processId, int machineId, int time) {
        result.machineWorkTime[machineId] = result.startTime[jobId][processId] + time;
    }

    // 更新作业的结束时间
    private void updateEndTime(Result result, int jobId, int processId, int machineId) {
        result.endTime[jobId][processId] = result.machineWorkTime[machineId]; // 更新作业的结束时间
        result.processIds[jobId] += 1; // 更新作业的处理过程ID
    }



    // 交叉基因
    private Gene crossGene(Gene g1, Gene g2) {
        List<Integer> indexList = makeList(chromosome_Size); // 生成索引列表
        int p1 = indexList.remove(random.nextInt(indexList.size())); // 随机选择交叉点1
        int p2 = indexList.remove(random.nextInt(indexList.size())); // 随机选择交叉点2
        int start = Math.min(p1, p2); // 计算交叉起始点
        int end = Math.min(Math.max(p1, p2), chromosome_Size - 1); // 计算交叉结束点，确保不超出染色体范围

        List<Integer> proto = subArray(g1.chromosome, start, end + 1); // 从第一个基因中获取交叉片段
        List<Integer> t = new ArrayList<>(); // 创建一个空列表用于存储第二个基因的染色体

        // 复制第二个基因的染色体
        for (Integer c : g2.chromosome) {
            t.add(c);
        }

        // 从第二个基因中移除交叉片段中的元素
        for (Integer val : proto) {
            for (int i = 0; i < t.size(); i++) {
                if (val.equals(t.get(i))) {
                    t.remove(i);
                    break;
                }
            }
        }

        Gene child = new Gene(); // 创建一个新的基因对象作为子代
        int sublistEndIndex = Math.min(chromosome_Size - proto.size(), t.size()); // 计算子列表的结束索引

        proto.addAll(t.subList(0, sublistEndIndex)); // 将交叉片段和剩余基因片段合并
        List<Integer> temp = new ArrayList<>(t.subList(sublistEndIndex, t.size()));
        temp.addAll(proto);
        child.chromosome = temp.toArray(new Integer[0]); // 将列表转换为数组作为子代的染色体
        child.fitness = (double) calculateFitness(child).fulfillTime; // 计算子代的适应度
        return child; // 返回交叉后的子代基因
    }



    // 变异基因
    public Gene mutationGene(Gene gene, int n) {
        // 生成染色体索引列表
        List<Integer> indexList = makeList(gene.chromosome.length);

        // 对基因进行n次变异
        for (int i = 0; i < n && !indexList.isEmpty(); i++) {
            int a = indexList.remove(random.nextInt(indexList.size())); // 随机选择变异位置a
            int b = indexList.isEmpty() ? a : indexList.remove(random.nextInt(indexList.size())); // 随机选择变异位置b
            int t = gene.chromosome[a];
            gene.chromosome[a] = gene.chromosome[b];
            gene.chromosome[b] = t;
        }

        // 重新计算基因的适应度
        gene.fitness = calculateFitness(gene).fulfillTime;

        return gene; // 返回变异后的基因
    }

    // 选择基因
    public Gene selectGene(int n) {
        // 生成基因集合的索引列表
        List<Integer> indexList = makeList(geneSet.size());

        // 随机选择n个索引作为被选中的基因
        Map<Integer, Boolean> map = new HashMap<>();
        for (int i = 0; i < n; i++) {
            map.put(indexList.remove(random.nextInt(indexList.size())), true);
        }

        // 初始化最佳基因为一个适应度极大值
        Gene bestGene = new Gene(0xfffff);
        int i = 0;

        // 遍历基因集合，找到被选中的最佳基因
        for (Gene gene : geneSet) {
            if (map.containsKey(i)) {
                if (bestGene.fitness > gene.fitness) {
                    bestGene = gene;
                }
            }
            i++;
        }

        return bestGene; // 返回被选中的最佳基因
    }


    // 运行遗传算法
    public Result run(List<List<int[]>> job) {
        int jobSize = job.size();

        // 检查工作数量是否超出范围
        if (jobSize > job_num) {
            System.out.println("工作数量超出范围！");
            return null;
        }

        // 处理作业信息
        for (int i = 0; i < jobSize; i++) {
            chromosome_Size += job.get(i).size();
            process_num = Math.max(process_num, job.get(i).size());
            for (int j = 0; j < job.get(i).size(); j++) {
                machine_Matrix[i][j] = job.get(i).get(j)[0];
                time_Matrix[i][j] = job.get(i).get(j)[1];
            }
        }

        // 构建工序矩阵
        for (int i = 0; i < jobSize; i++) {
            for (int j = 0; j < process_num; j++) {
                if (machine_Matrix[i][j] != -1) {
                    process_Matrix[i][machine_Matrix[i][j]] = j;
                }
            }
        }

        // 初始化种群
        initialPopulation();

        // 运行遗传算法
        for (int i = 0; i < population_num; i++) {
            double p = (double) random.nextInt(100) / 100.0;
            if (p < rate) {
                int index = random.nextInt(geneSet.size());
                int k = 0;
                for (Gene gene : geneSet) {
                    if (k == index) {
                        mutationGene(gene, 2); // 变异基因
                        break;
                    }
                    k++;
                }
            } else {
                Gene g1 = selectGene(3), g2 = selectGene(3);
                Gene child1 = crossGene(g1, g2); // 交叉基因
                Gene child2 = crossGene(g2, g1);
                geneSet.add(child1);
                geneSet.add(child2);
            }
        }

        // 找到适应度最佳的基因
        Gene bestGene = new Gene(0xffffff);
        for (Gene gene : geneSet) {
            if (bestGene.fitness > gene.fitness) {
                bestGene = gene;
            }
        }

        return calculateFitness(bestGene); // 返回最优解的适应度
    }

    // 选择基因（默认选择3个）
    public Gene selectGene() {
        return selectGene(3);
    }

    // 变异基因（默认变异2次）
    public Gene mutationGene(Gene gene) {
        return mutationGene(gene, 2);
    }

    // 主函数
    public static void main(String[] args) {

        int n = 4; // 机器数
        int m = 5; // 工件数

        List<List<int[]>> job = new ArrayList<>();
        job.add(Arrays.asList(new int[]{0, 3}, new int[]{1, 2}, new int[]{2, 2}, new int[]{3, 4}, new int[]{1, 3}));
        job.add(Arrays.asList(new int[]{0, 2}, new int[]{2, 1}, new int[]{1, 4}, new int[]{3, 2}, new int[]{0, 3}));
        job.add(Arrays.asList(new int[]{1, 4}, new int[]{2, 3}, new int[]{0, 2}, new int[]{3, 3}, new int[]{1, 2}));
        job.add(Arrays.asList(new int[]{2, 3}, new int[]{1, 2}, new int[]{3, 1}, new int[]{0, 4}, new int[]{2, 3}));

        GeneticAlgorithm ga = new GeneticAlgorithm(n, m);
        Result result = ga.run(job);
        int processNumber = ga.process_num;
        int[][] machineMatrix = ga.machine_Matrix;

        System.out.println(String.format("最优解: %d",result.fulfillTime));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < processNumber; j++) {
                if (machineMatrix[i][j] != -1) {
                    System.out.println(String.format("job: %d, process: %d, machine: %d, startTime: %d, endTime: %d",
                            i, j, machineMatrix[i][j], result.startTime[i][j], result.endTime[i][j]));
                }
            }
        }
    }
}