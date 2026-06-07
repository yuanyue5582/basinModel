#include "fillDEM.h"
#include "matrixData.h"
#include "utils.h"
#include "HashPriorityQueue.h"
#include <iostream>
#include <chrono>
#include <cmath>
#include <queue>
#include <vector>
#include <cstring>
#include "gdal_priv.h"

// ============================================================
// 辅助函数（修正：设置流向前检查是否已存在）
// ============================================================
static inline void setFlowDirection(std::vector<uint8_t>& flow_direction,
                                    int from_row, int from_col,
                                    int to_row, int to_col,
                                    int width) {
    int from_idx = from_row * width + from_col;
    // 如果已有流向，不再覆盖
    if (flow_direction[from_idx] != 0) return;

    int dr = to_row - from_row;
    int dc = to_col - from_col;
    int dir_idx = -1;
    if (dr == 0 && dc == 1) dir_idx = 0;
    else if (dr == 1 && dc == 1) dir_idx = 1;
    else if (dr == 1 && dc == 0) dir_idx = 2;
    else if (dr == 1 && dc == -1) dir_idx = 3;
    else if (dr == 0 && dc == -1) dir_idx = 4;
    else if (dr == -1 && dc == -1) dir_idx = 5;
    else if (dr == -1 && dc == 0) dir_idx = 6;
    else if (dr == -1 && dc == 1) dir_idx = 7;
    if (dir_idx != -1) {
        static constexpr uint8_t dir_code[8] = {1, 2, 4, 8, 16, 32, 64, 128};
        flow_direction[from_idx] = dir_code[dir_idx];
    }
}

// 标记集水区，但仅当邻居无流向时才进行（新要求）
static inline void markCatchmentIfNeeded(std::vector<uint32_t>& catchment_matrix,
                                         int idx, uint32_t current_label,
                                         const std::vector<uint8_t>& flow_direction) {
    if (current_label != 0 && catchment_matrix[idx] == 0 && flow_direction[idx] == 0) {
        catchment_matrix[idx] = current_label;
    }
}

// ============================================================
// 主处理流程
// ============================================================
bool FillDEM::process(const std::string& input_dem, const std::string& output_dir) {
    std::cout << "=== Module 1: Barnes Fill (No DEM modification) + Tree Building ===" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    // 1. 读取DEM
    GDALAllRegister();
    GDALDataset* poDataset = (GDALDataset*)GDALOpen(input_dem.c_str(), GA_ReadOnly);
    if (!poDataset) {
        std::cerr << "无法打开DEM文件: " << input_dem << std::endl;
        return false;
    }
    width = poDataset->GetRasterXSize();
    height = poDataset->GetRasterYSize();
    double geo_transform[6];
    poDataset->GetGeoTransform(geo_transform);
    cell_width = std::abs(static_cast<float>(geo_transform[1]));
    cell_height = std::abs(static_cast<float>(geo_transform[5]));

    std::vector<float> dem_data(width * height);
    GDALRasterBand* poBand = poDataset->GetRasterBand(1);
    int pbSuccess = 0;
    double dfNoData = poBand->GetNoDataValue(&pbSuccess);
    if (pbSuccess) no_data_value = static_cast<float>(dfNoData);
    poBand->RasterIO(GF_Read, 0, 0, width, height,
                     dem_data.data(), width, height, GDT_Float32, 0, 0);
    GDALClose((GDALDatasetH)poDataset);
    std::cout << "DEM尺寸: " << width << " x " << height << std::endl;

    // 2. 初始化矩阵
    std::vector<uint32_t> depression_matrix(width * height, 0);
    std::vector<uint8_t> flow_direction(width * height, 0);
    std::vector<uint32_t> catchment_matrix(width * height, 0);

    // 3. 执行填洼算法（仅标记洼地，不修改DEM）
    auto phase1_start = std::chrono::high_resolution_clock::now();
    std::cout << "执行Barnes填洼算法（只记录洼地，不填平）..." << std::endl;
    fillDEM_Barnes(dem_data, depression_matrix, flow_direction, catchment_matrix);
    auto phase1_end = std::chrono::high_resolution_clock::now();
    auto phase1_ms = std::chrono::duration_cast<std::chrono::milliseconds>(phase1_end - phase1_start).count();
    std::cout << "第一阶段耗时: " << phase1_ms << " ms" << std::endl;

    // 4. 直接从 pending_outlets 构建树节点并建立父子关系
    auto build_start = std::chrono::high_resolution_clock::now();
    std::cout << "从洼地出口信息直接构建树节点及关系..." << std::endl;

    int32_t max_dep_id = 0;
    for (const auto& info : pending_outlets) {
        if (info.dep_id > max_dep_id) max_dep_id = info.dep_id;
    }
    tree.reserveNodes(max_dep_id);
    for (const auto& info : pending_outlets) {
        tree.addNode(info.dep_id, info.row, info.col);
        if (info.parent != -1) {
            tree.setParent(info.dep_id, info.parent);
        } else {
            tree.markBoundary(info.dep_id);
        }
    }
    tree.finalize();

    auto build_end = std::chrono::high_resolution_clock::now();
    auto build_ms = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start).count();
    std::cout << "构建树节点及关系耗时: " << build_ms << " ms" << std::endl;

    // 5. 保存结果（不保存填洼DEM）
    if (!createDirectory(output_dir)) {
        std::cerr << "无法创建输出目录: " << output_dir << std::endl;
        return false;
    }
    MatrixData mat;
    mat.setDimensions(width, height);
    mat.setGeoTransform(geo_transform);

    std::string dep_id_path = output_dir + "/depression_id.tif";
    std::string flow_path = output_dir + "/flow_direction.tif";
    std::string catch_path = output_dir + "/catchment_areas.tif";
    std::string tree_path = output_dir + "/hash_weight_tree.bin";

    if (!mat.saveDepressionID(dep_id_path, depression_matrix)) {
        std::cerr << "保存洼地ID矩阵失败!" << std::endl;
        return false;
    }
    if (!mat.saveFlowDirection(flow_path, flow_direction)) {
        std::cerr << "保存流向矩阵失败!" << std::endl;
        return false;
    }
    if (!mat.saveCatchmentID(catch_path, catchment_matrix)) {
        std::cerr << "保存集水区矩阵失败!" << std::endl;
        return false;
    }
    std::cout << "洼地ID矩阵已保存: " << dep_id_path << std::endl;
    std::cout << "流向矩阵已保存: " << flow_path << std::endl;
    std::cout << "集水区矩阵已保存: " << catch_path << std::endl;

    if (!tree.saveToBinaryFile(tree_path)) {
        std::cerr << "保存哈希权重树失败!" << std::endl;
        return false;
    }
    std::cout << "哈希权重树已保存: " << tree_path << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "\n模块1完成！" << std::endl;
    std::cout << "洼地数量: " << tree.getNodeCount() << std::endl;
    std::cout << "总耗时: " << duration.count() << " ms" << std::endl;

    return true;
}

// ============================================================
// Barnes填洼算法（修正：集水区传播仅当邻居无流向）
// ============================================================
void FillDEM::fillDEM_Barnes(const std::vector<float>& dem_data,
                             std::vector<uint32_t>& depression_matrix,
                             std::vector<uint8_t>& flow_direction,
                             std::vector<uint32_t>& catchment_matrix) {
    HashPriorityQueue priority_queue;
    std::queue<Node> plain_queue;
    // processed[i] = 1 表示该单元格已经处理完毕（邻居已处理，流向已确定）
    std::vector<uint8_t> processed(width * height, 0);

    // 初始化：将所有边界单元格加入优先队列（不标记为已处理）
    std::vector<std::pair<int, int>> boundaries;
    findValidBoundaries(dem_data, boundaries);
    for (const auto& boundary : boundaries) {
        int row = boundary.first, col = boundary.second;
        if (!processed[row * width + col]) {
            priority_queue.push(Node(row, col, dem_data[row * width + col]));
        }
    }
    std::cout << "边界单元格数: " << boundaries.size() << std::endl;

    int processed_cells = 0;
    int created_depressions = 0;
    int next_dep_id = 1;
    pending_outlets.clear();

    // 外层循环：优先队列弹出
    while (!priority_queue.empty()) {
        Node current_node = priority_queue.pop();
        int c_row = current_node.row;
        int c_col = current_node.col;
        float c_value = current_node.spill;   // 原始高程
        int c_idx = c_row * width + c_col;

        // 如果该单元格已经处理过，则跳过（避免重复处理）
        if (processed[c_idx]) continue;

        int current_dep_id = depression_matrix[c_idx];

        // 遍历八个邻居
        for (int i = 0; i < 8; i++) {
            int n_row = Get_rowTo(i, c_row);
            int n_col = Get_colTo(i, c_col);
            if (!is_InGrid(n_row, n_col)) continue;
            int n_idx = n_row * width + n_col;
            if (processed[n_idx]) continue;   // 只跳过已处理的邻居

            float neighbor_value = dem_data[n_idx];
            if (is_NoData(neighbor_value)) continue;

            // 邻居更低 -> 洼地扩展或创建新洼地
            if (neighbor_value < c_value) {
                // 若当前单元格还不是洼地，则新建洼地（当前单元格作为出口）
                if (current_dep_id == 0) {
                    current_dep_id = next_dep_id++;
                    created_depressions++;
                    depression_matrix[c_idx] = current_dep_id;
                    catchment_matrix[c_idx] = current_dep_id;

                    // 获取出口流向：优先使用已有，否则计算D8
                    uint8_t outlet_dir = flow_direction[c_idx];
                    if (outlet_dir == 0) {
			//std::cout << current_dep_id << " 的出口 " << std::endl;
                        outlet_dir = computeD8Flow(c_row, c_col, dem_data);
                        if (outlet_dir != 0 && flow_direction[c_idx] == 0) {
                            flow_direction[c_idx] = outlet_dir;
                        }
                    }

                    // 确定父洼地ID（下游单元格的洼地/集水区ID）
                    int32_t parent_id = -1;
                    if (outlet_dir != 0) {
                        int dir_idx = -1;
                        for (int k = 0; k < 8; ++k) {
                            if (dir_code[k] == outlet_dir) {
                                dir_idx = k;
                                break;
                            }
                        }
                        if (dir_idx != -1) {
                            int down_row = Get_rowTo(dir_idx, c_row);
                            int down_col = Get_colTo(dir_idx, c_col);
                            if (is_InGrid(down_row, down_col)) {
                                int down_idx = down_row * width + down_col;
                                uint32_t target = depression_matrix[down_idx];
                                if (target == 0) target = catchment_matrix[down_idx];
                                if (target != 0) parent_id = static_cast<int32_t>(target);
                            }
                        }
                    }
                    pending_outlets.push_back({current_dep_id, c_row, c_col, parent_id});

                    // 检查周围是否有流向指向当前单元格的邻居，将其集水区标记为当前洼地
                    for (int j = 0; j < 8; j++) {
                        int nb_row = Get_rowTo(j, c_row);
                        int nb_col = Get_colTo(j, c_col);
                        if (!is_InGrid(nb_row, nb_col)) continue;
                        int nb_idx = nb_row * width + nb_col;
                        if (depression_matrix[nb_idx] == 0 && catchment_matrix[nb_idx] == 0 && flow_direction[nb_idx] != 0) {
                            int target_row = nb_row, target_col = nb_col;
                            bool found = false;
                            for (int k = 0; k < 8; k++) {
                                if (dir_code[k] == flow_direction[nb_idx]) {
                                    target_row = Get_rowTo(k, nb_row);
                                    target_col = Get_colTo(k, nb_col);
                                    found = true;
                                    break;
                                }
                            }
                            if (found && target_row == c_row && target_col == c_col) {
                                catchment_matrix[nb_idx] = current_dep_id;
                            }
                        }
                    }
                }

                // 将邻居加入普通队列（洼地内部单元格），不标记 processed
                if (depression_matrix[n_idx] == 0) {
                    depression_matrix[n_idx] = current_dep_id;
                    catchment_matrix[n_idx] = current_dep_id;
                    plain_queue.push(Node(n_row, n_col, c_value));
                }
            }
            // 邻居更高或相等 -> 加入优先队列，并设置流向（仅当邻居无流向时）
            else {
                // 只有当邻居还没有流向时，才设置其流向并标记集水区
                if (flow_direction[n_idx] == 0) {
                    setFlowDirection(flow_direction, n_row, n_col, c_row, c_col, width);
                    uint32_t current_label = (current_dep_id != 0) ? current_dep_id : catchment_matrix[c_idx];
                    if (current_label != 0 && catchment_matrix[n_idx] == 0) {
                        catchment_matrix[n_idx] = current_label;
                    }
                   priority_queue.push(Node(n_row, n_col, neighbor_value));
                }
            }
        }

        // 当前单元格处理完毕，标记为已处理
        processed[c_idx] = 1;
        processed_cells++;

        // ========== 立即处理普通队列中所有洼地内部单元格 ==========
        while (!plain_queue.empty()) {
            Node plain_node = plain_queue.front();
            plain_queue.pop();
            int p_row = plain_node.row;
            int p_col = plain_node.col;
            float p_spill = plain_node.spill;   // 出口高程（阈值）
            int p_idx = p_row * width + p_col;

            // 如果该单元格已经处理过，则跳过
            if (processed[p_idx]) continue;

            int p_dep_id = depression_matrix[p_idx];   // 该洼地的ID
            processed_cells++;

            // 遍历邻居
            for (int i = 0; i < 8; i++) {
                int n_row = Get_rowTo(i, p_row);
                int n_col = Get_colTo(i, p_col);
                if (!is_InGrid(n_row, n_col)) continue;
                int n_idx = n_row * width + n_col;
                if (processed[n_idx]) continue;

                float neighbor_value = dem_data[n_idx];
                if (is_NoData(neighbor_value)) continue;

                if (neighbor_value < p_spill) {
                    // 邻居更低，继续扩展同一洼地
                    if (depression_matrix[n_idx] == 0) {
                        depression_matrix[n_idx] = p_dep_id;
                        catchment_matrix[n_idx] = p_dep_id;
                        plain_queue.push(Node(n_row, n_col, p_spill));
                    }
                } else {
                    // 邻居更高或相等，加入优先队列（仅当邻居无流向时设置流向和集水区）
                    if (flow_direction[n_idx] == 0) {
                        setFlowDirection(flow_direction, n_row, n_col, p_row, p_col, width);
                        if (p_dep_id != 0 && catchment_matrix[n_idx] == 0) {
                            catchment_matrix[n_idx] = p_dep_id;
                        }
                        priority_queue.push(Node(n_row, n_col, neighbor_value));
                    }
                }
            }

            // 当前洼地内部单元格处理完毕，标记为已处理
            processed[p_idx] = 1;
        } // while plain_queue
    } // while priority_queue

    std::cout << "\n处理完成！" << std::endl;
    std::cout << "总处理单元格: " << processed_cells << std::endl;
    std::cout << "创建洼地数: " << created_depressions << std::endl;
}

// ============================================================
// 辅助函数：找到所有边界单元格（用于初始化优先队列）
// ============================================================
bool FillDEM::findValidBoundaries(const std::vector<float>& dem_data,
                                  std::vector<std::pair<int,int>>& boundaries) {
    boundaries.clear();
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            if (is_NoData(dem_data[row * width + col])) continue;
            bool isBoundary = false;
            for (int i = 0; i < 8; i++) {
                int nRow = Get_rowTo(i, row);
                int nCol = Get_colTo(i, col);
                if (!is_InGrid(nRow, nCol) || is_NoData(dem_data[nRow * width + nCol])) {
                    isBoundary = true;
                    break;
                }
            }
            if (isBoundary) boundaries.emplace_back(row, col);
        }
    }
    return !boundaries.empty();
}

// ============================================================
// D8 最陡坡降流向计算
// ============================================================
uint8_t FillDEM::computeD8Flow(int row, int col, const std::vector<float>& dem_data) const {
    float cur_elev = dem_data[row * width + col];
    float max_slope = 0.0f;
    uint8_t best_dir = 0;
    for (int i = 0; i < 8; i++) {
        int n_row = Get_rowTo(i, row);
        int n_col = Get_colTo(i, col);
        if (!is_InGrid(n_row, n_col)) continue;
        float n_elev = dem_data[n_row * width + n_col];
        if (is_NoData(n_elev)) continue;
        if (n_elev >= cur_elev) continue;
        float drop = cur_elev - n_elev;
        float distance = (i % 2 == 0) ? 1.0f : 1.41421356f;
        float slope = drop / distance;
        if (slope > max_slope) {
            max_slope = slope;
            best_dir = dir_code[i];
        }
    }
    return best_dir;
}

// ============================================================
// 获取从 (from_r, from_c) 到 (to_r, to_c) 的方向索引（0-7）
// ============================================================
int FillDEM::getDirectionIndex(int from_r, int from_c, int to_r, int to_c) const {
    int dr = to_r - from_r;
    int dc = to_c - from_c;
    if (dr == 0 && dc == 1) return 0;
    if (dr == 1 && dc == 1) return 1;
    if (dr == 1 && dc == 0) return 2;
    if (dr == 1 && dc == -1) return 3;
    if (dr == 0 && dc == -1) return 4;
    if (dr == -1 && dc == -1) return 5;
    if (dr == -1 && dc == 0) return 6;
    if (dr == -1 && dc == 1) return 7;
    return -1;
}