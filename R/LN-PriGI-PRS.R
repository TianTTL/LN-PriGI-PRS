#!/usr/bin/env Rscript

# LN-PriGI-PRS (Lupus nephritis specific Prioritized-Gene-Informed PRS)
# data: 2026-05-29
# author: Xingjian Gao

# 加载必要的R包
suppressPackageStartupMessages({
    library(optparse)
    library(this.path)
    library(data.table)
    library(dplyr)
})

# 加载自定义模块
source(file.path(this.dir(), "preprocess.R"))
source(file.path(this.dir(), "run_plink_pgs.R"))

# 解析命令行参数
# @return 解析后的参数列表
parse_arguments <- function() {
    option_list <- list(
        make_option(
            c("-i", "--input"), 
            type = "character", 
            default = NULL,
            help = "用户输入的PLINK文件前缀（必需）"),
        
        make_option(
            c("-o", "--output"), 
            type = "character", 
            default = NULL,
            help = "结果输出文件路径（必需）"),
        
        make_option(
            c("--overlap_snp_cutoff"), 
            type = "double", 
            default = 0.9,
            help = "SNP重叠率阈值 [默认: %default]"),
        
        make_option(
            c("--data_dir"), 
            type = "character", 
            default = NULL,
            help = "内置数据目录路径 [默认: 自动检测]"),
        
        make_option(
            c("--temp_dir"), 
            type = "character", 
            default = NULL,
            help = "临时文件目录 [默认: 输出目录]"),
        
        make_option(
            c("--version"), 
            action = "store_true", 
            default = FALSE,
            help = "显示版本信息")
    )
    
    parser <- OptionParser(
        option_list = option_list, 
        description = "狼疮性肾炎多基因评分计算系统")
    
    args <- parse_args(parser)
    
    # 检查必需参数
    if (is.null(args$input)) {
        cat("错误: 必须指定输入文件前缀 (--input)\n")
        print_help(parser)
        quit(status = 1)
    }
    
    if (is.null(args$output)) {
        cat("错误: 必须指定输出文件路径 (--output)\n")
        print_help(parser)
        quit(status = 1)
    }
    
    # 显示版本信息
    if (args$version) {
        cat("狼疮性肾炎多基因评分计算系统 v1.0\n")
        cat("基于R语言和plink的狼疮性肾炎多基因评分计算系统\n")
        quit(status = 0)
    }
    
    return(args)
}

# 检查系统依赖
# @return 无返回值，如果检查失败则退出程序
check_dependencies <- function() {
    cat("检查系统依赖...\n")
    
    # 检查PLINK是否可用
    plink_check <- system("plink --version", ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (plink_check != 0) {
        cat("错误: PLINK未安装或不在PATH中\n")
        cat("请安装PLINK并确保其在系统PATH中\n")
        cat("下载地址: https://www.cog-genomics.org/plink2/\n")
        quit(status = 1)
    }
    
    # 检查R包是否可用
    required_packages <- c("optparse", "this.path", "data.table", "dplyr")
    missing_packages <- c()
    
    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            missing_packages <- c(missing_packages, pkg)
        }
    }
    
    if (length(missing_packages) > 0) {
        cat("错误: 缺少必需的R包:", paste(missing_packages, collapse = ", "), "\n")
        cat("请使用以下命令安装:\n")
        cat("install.packages(c(", paste0("'", missing_packages, "'", collapse = ", "), "))\n")
        quit(status = 1)
    }
    
    cat("系统依赖检查通过\n")
}

# 清理临时文件
# @param temp_files 临时文件列表
cleanup_temp_files <- function(temp_files) {
    for (file in temp_files) {
        if (file.exists(file)) {
            file.remove(file)
        }
    }
}

# 主函数
main <- function() {
    # 解析命令行参数
    args <- parse_arguments()

    ## 检查SNP重叠率阈值
    if (args$overlap_snp_cutoff < 0 || args$overlap_snp_cutoff > 1) {
        stop("SNP重叠率阈值必须在0到1之间")
    }
    
    # 显示程序信息
    cat("========================================\n")
    cat("狼疮性肾炎多基因评分计算系统 v1.0.0\n")
    cat("========================================\n")
    
    # 检查系统依赖
    check_dependencies()
    
    # 设置数据目录
    if (is.null(args$data_dir)) {
        args$data_dir <- file.path(dirname(this.dir()), "data")
    }
    
    # 设置临时文件目录
    if (is.null(args$temp_dir)) {
        args$temp_dir <- dirname(args$output)
    }
    
    # 创建输出目录
    output_dir <- dirname(args$output)
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # 创建临时文件目录
    if (!dir.exists(args$temp_dir)) {
        dir.create(args$temp_dir, recursive = TRUE)
    }
    
    # 临时文件列表
    temp_files <- c()
    
    tryCatch({
        # 步骤1: 数据预处理
        cat("\n步骤1: 数据预处理\n")
        cat("==================\n")
              
        preprocess_result <- preprocess_data(
            input_prefix = args$input,
            output_prefix = args$temp_dir,
            data_dir = args$data_dir
        )
        
        # 步骤2: PGS计算
        cat("\n步骤2: PGS计算\n")
        cat("===============\n")

        pgs_temp_prefix <- file.path(args$temp_dir, "temp_pgs")
        temp_files <- c(
            temp_files, 
            paste0(pgs_temp_prefix, c(".profile", ".log", ".nosex")),
            file.path(args$temp_dir, "temp_weights.txt"))

        pgs_results <- calculate_pgs(
            input_prefix = args$input,
            matched_snps = preprocess_result$matched_snps,
            output_prefix = pgs_temp_prefix,
            temp_dir = args$temp_dir
        )
        
        # 步骤3: 分位数计算
        cat("\n步骤3: 分位数计算\n")
        cat("==================\n")
        
        pgs_with_percentiles <- calc_percentile(
            pgs_results = pgs_results,
            reference_pgs = preprocess_result$model_data$reference_pgs
        )
        
        # 步骤4: 生成最终报告
        cat("\n步骤4: 生成最终报告\n")
        cat("====================\n")
        
        final_results <- generate_pgs_report(
            pgs_results = pgs_with_percentiles,
            snp_overlap_rate = preprocess_result$snp_overlap_rate,
            overlap_snp_cutoff = args$overlap_snp_cutoff
        )
        
        # 重新排列列顺序
        final_results <- final_results %>%
            select(SAMPLE_ID, PGS_SCORE, PERCENTILE, RISK_LEVEL, WARNING)
        
        # 写入结果文件
        fwrite(final_results, paste0(args$output, ".report"), sep = "\t", quote = FALSE)
        
        # 显示结果摘要
        cat("\n结果摘要\n")
        cat("========\n")
        cat("输出文件:", args$output, "\n")
        cat("样本数量:", nrow(final_results), "\n")
        cat("使用SNP数:", preprocess_result$used_snps, "/", preprocess_result$total_model_snps, "\n")
        cat("SNP重合率:", round(preprocess_result$snp_overlap_rate * 100, 2), "%\n")
        
        if (preprocess_result$snp_overlap_rate < args$overlap_snp_cutoff) {
            cat("警告: SNP重合率过低，结果可能不可靠\n")
        }
        
        cat("\nPGS得分范围:", round(range(final_results$PGS_SCORE, na.rm = TRUE), 3), "\n")
        cat("分位数范围:", round(range(final_results$PERCENTILE, na.rm = TRUE), 3), "\n")
        
        cat("\n分析完成！\n")
        
    }, error = function(e) {
        cat("错误:", e$message, "\n")
        quit(status = 1)
    }, finally = {
        # 清理临时文件
        if (length(temp_files) > 0) {
            cat("\n清理临时文件...\n")
            cleanup_temp_files(temp_files)
        }
    })
}

# 如果直接运行此脚本，则执行主函数
if (!interactive()) {
    main()
}

# 狼疮性肾炎多基因评分计算系统 - 数据预处理模块
# 功能：数据加载、校验、预处理、质控

# 加载必要的R包
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

# 获取脚本所在目录的绝对路径
# @return 脚本所在目录的绝对路径
get_script_dir <- function() {
    # 获取当前脚本的路径
    script_path <- commandArgs(trailingOnly = FALSE)
    script_path <- script_path[grepl("--file=", script_path)]
    if (length(script_path) > 0) {
        script_path <- sub("--file=", "", script_path)
        return(dirname(normalizePath(script_path)))
    } else {
        # 如果无法获取脚本路径，使用当前工作目录
        return(getwd())
    }
}

# 加载内置模型数据
# @param data_dir 数据目录路径
# @return 包含模型SNP、权重和参考数据的列表
load_model_data <- function(data_dir) {
    cat("正在加载内置模型数据...\n")
    
    # 检查数据目录是否存在
    if (!dir.exists(data_dir)) {
        stop("数据目录不存在: ", data_dir)
    }
    
    # 加载模型SNP列表
    model_snps_file <- file.path(data_dir, "model_snps.txt")
    if (!file.exists(model_snps_file)) {
        stop("模型SNP文件不存在: ", model_snps_file)
    }
    model_snps <- fread(model_snps_file, header = FALSE, skip = 1)
    colnames(model_snps) <- c("SNP_ID", "CHR", "POS", "A1", "A2", "A1Effect")
    cat("加载了", nrow(model_snps), "个模型SNP\n")
    
    # 加载参考人群PGS
    reference_file <- file.path(data_dir, "reference_pgs.txt")
    if (!file.exists(reference_file)) {
        stop("参考人群PGS文件不存在: ", reference_file)
    }
    reference_pgs <- fread(reference_file, header = FALSE)
    colnames(reference_pgs) <- c("SAMPLE_ID", "PGS_SCORE")
    cat("加载了", nrow(reference_pgs), "个参考样本的PGS数据\n")
    
    return(list(
        model_snps = model_snps,
        reference_pgs = reference_pgs
    ))
}

# 加载用户PLINK数据
# @param input_prefix PLINK文件前缀
# @return 包含.bim和.fam数据的列表
load_user_data <- function(input_prefix) {
    cat("正在加载用户PLINK数据...\n")
    
    # 检查PLINK文件是否存在
    bim_file <- paste0(input_prefix, ".bim")
    fam_file <- paste0(input_prefix, ".fam")
    bed_file <- paste0(input_prefix, ".bed")
    
    if (!file.exists(bim_file)) {
        stop("PLINK .bim文件不存在: ", bim_file)
    }
    if (!file.exists(fam_file)) {
        stop("PLINK .fam文件不存在: ", fam_file)
    }
    if (!file.exists(bed_file)) {
        stop("PLINK .bed文件不存在: ", bed_file)
    }
    
    # 加载.bim文件
    bim_data <- fread(bim_file, header = FALSE)
    colnames(bim_data) <- c("CHR", "SNP_ID", "GENETIC_DIST", "POS", "ALT_ALLELE", "REF_ALLELE")
    
    # 加载.fam文件
    fam_data <- fread(fam_file, header = FALSE)
    colnames(fam_data) <- c("FAM_ID", "IND_ID", "PAT_ID", "MAT_ID", "SEX", "PHENOTYPE")

    cat("用户数据包含", nrow(bim_data), "个SNP和", nrow(fam_data), "个样本\n")

    return(list(
        bim_data = bim_data,
        fam_data = fam_data
    ))
}

# SNP匹配与等位基因对齐
# @param user_bim 用户.bim数据
# @param model_snps 模型SNP数据
# @return 对齐后的SNP列表
align_snps <- function(user_bim, model_snps) {
    cat("正在进行SNP匹配与等位基因对齐...\n")
    
    # 根据染色体和位置匹配SNP
    user_bim$CHR_POS <- paste(user_bim$CHR, user_bim$POS, sep = "_")
    model_snps$CHR_POS <- paste(model_snps$CHR, model_snps$POS, sep = "_")
    
    # 找到匹配的SNP
    matched_snps <- left_join(user_bim, model_snps, by = "CHR_POS") %>%
        filter(!is.na(SNP_ID.y)) %>%
        filter(ALT_ALLELE == A1 | REF_ALLELE == A1) %>%
        select(SNP_ID.x, A1, A2, A1Effect) %>%
        rename(SNP_ID = SNP_ID.x)

    cat("找到", nrow(matched_snps), "个匹配的SNP\n")
    
    if (nrow(matched_snps) == 0) {
        stop("没有找到匹配的SNP，请检查数据格式和坐标系统")
    }

    return(matched_snps)
}

# 获取等位基因的互补碱基
# @param allele 等位基因
# @return 互补等位基因
complement_allele <- function(allele) {
    complement_map <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G")
    return(complement_map[allele])
}

# 主预处理函数
# @param input_prefix 用户PLINK文件前缀
# @param output_prefix 输出文件前缀
# @param data_dir 内置数据目录
# @return 包含质控后文件前缀和SNP信息的列表
preprocess_data <- function(input_prefix, output_prefix, data_dir = NULL) {
    
    # 加载模型数据
    model_data <- load_model_data(data_dir)
    
    # 加载用户数据
    user_data <- load_user_data(input_prefix)
    
    # SNP匹配与对齐
    matched_snps <- align_snps(user_data$bim_data, model_data$model_snps)
    
    # 计算SNP缺失率
    total_model_snps_n <- nrow(model_data$model_snps)
    used_snps_n <- nrow(matched_snps)
    snp_overlap_rate <- used_snps_n / total_model_snps_n
    
    cat("SNP使用情况统计:\n")
    cat("  模型总SNP数:", total_model_snps_n, "\n")
    cat("  实际使用SNP数:", used_snps_n, "\n")
    cat("  SNP重合率:", round(snp_overlap_rate * 100, 2), "%\n")
    
    return(list(
        matched_snps = matched_snps,
        snp_overlap_rate = snp_overlap_rate,
        total_model_snps_n = total_model_snps_n,
        used_snps_n = used_snps_n,
        model_data = model_data
    ))
}

# 狼疮性肾炎多基因评分计算系统 - PLINK PGS计算模块
# 功能：调用PLINK进行PGS计算

# 加载必要的R包
suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

# 创建PLINK权重文件
# @param matched_snps 对齐后的SNP数据
# @param temp_weights_file 输出权重文件路径
create_plink_weights_file <- function(matched_snps, temp_weights_file) {
    cat("正在创建PLINK权重文件...\n")
    
    # 创建PLINK格式的权重文件
    # 格式：SNP_ID 等位基因 权重
    plink_weights <- data.table(
        SNP_ID = matched_snps$SNP_ID,
        ALLELE = matched_snps$A1,  # PLINK使用ALT等位基因作为效应等位基因
        WEIGHT = matched_snps$A1Effect
    )
    
    # 写入文件
    fwrite(
        plink_weights, 
        temp_weights_file, 
        sep = "\t", col.names = TRUE, quote = FALSE)
    
    cat("权重文件已创建:", temp_weights_file, "\n")
    cat("包含", nrow(plink_weights), "个SNP的权重信息\n")
}

# 执行PLINK PGS计算
# @param input_prefix 用户输入的PLINK文件前缀
# @param temp_weights_file 权重文件路径
# @param output_prefix 输出文件前缀
# @return PGS结果文件路径
run_plink_score <- function(input_prefix, temp_weights_file, output_prefix) {
    cat("正在执行PLINK PGS计算...\n")
    
    # 构建PLINK score命令
    plink_cmd <- sprintf(
        "plink --bfile %s --score %s 1 2 3 header sum --out %s",
        input_prefix, temp_weights_file, output_prefix
    )
    
    cat("执行PLINK命令:", plink_cmd, "\n")
    
    # 执行PLINK命令
    exit_code <- system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if (exit_code != 0) {
        stop("PLINK PGS计算失败，退出代码:", exit_code)
    }
    
    # 检查输出文件是否存在
    score_file <- paste0(output_prefix, ".profile")
    if (!file.exists(score_file)) {
        stop("PLINK PGS计算结果文件不存在:", score_file)
    }
    
    cat("PLINK PGS计算完成\n")
    
    return(score_file)
}

# 解析PLINK PGS结果
# @param score_file PLINK score结果文件路径
# @return 包含样本ID和PGS得分的data.table
parse_plink_score <- function(score_file) {
    cat("正在解析PLINK PGS结果...\n")
    
    # 读取PLINK score结果
    # plink的.sscore文件格式：
    # #FID IID NMISS_ALLELE_CT DENOM SCORE1_AVG SCORE1_SUM
    score_data <- fread(score_file, header = TRUE)
    
    # 提取样本ID和PGS得分
    pgs_results <- score_data %>%
        select(IID, SCORESUM) %>%
        rename(SAMPLE_ID = IID, PGS_SCORE = SCORESUM)
    
    cat("解析了", nrow(pgs_results), "个样本的PGS结果\n")
    
    return(pgs_results)
}

# 主PGS计算函数
# @param input_prefix 用户输入的PLINK文件前缀
# @param matched_snps 对齐后的SNP数据
# @param output_prefix 输出文件前缀
# @param temp_dir 临时文件目录
# @return 包含PGS结果的data.table
calculate_pgs <- function(input_prefix, matched_snps, output_prefix, temp_dir = NULL) {
    # 创建临时权重文件
    temp_weights_file <- file.path(temp_dir, "temp_weights.txt")
    create_plink_weights_file(matched_snps, temp_weights_file)
    
    # 执行PLINK PGS计算
    score_file <- run_plink_score(input_prefix, temp_weights_file, output_prefix)
    
    # 解析PGS结果
    pgs_results <- parse_plink_score(score_file)
    
    return(pgs_results)
}

# 批量计算PGS分位数
# @param pgs_results 用户PGS结果
# @param reference_pgs 参考人群PGS数据
# @return 包含分位数的PGS结果
calc_percentile <- function(pgs_results, reference_pgs) {
    cat("正在计算PGS分位数...\n")
    
    # 为每个样本计算分位数
    pgs_results$PERCENTILE <- sapply(pgs_results$PGS_SCORE, function(score) {
        ecdf(reference_pgs$PGS_SCORE)(score)
    })
    
    cat("分位数计算完成\n")
    
    return(pgs_results)
}

# 生成PGS报告
# @param pgs_results 包含分位数的PGS结果
# @param snp_missing_rate SNP缺失率
# @param missing_snp_cutoff SNP缺失率阈值
# @return 处理后的PGS结果
generate_pgs_report <- function(pgs_results, snp_overlap_rate, overlap_snp_cutoff) {
    cat("正在生成PGS报告...\n")
    
    # 检查SNP缺失率
    if (snp_overlap_rate < overlap_snp_cutoff) {
        warning("SNP重合率 (", round(snp_overlap_rate * 100, 2), 
                        "%) 低于阈值 (", round(overlap_snp_cutoff * 100, 2), 
                        "%)，分位数结果可能不可靠")
        
        # 将分位数设为NA
        pgs_results$PERCENTILE <- NA
        pgs_results$WARNING <- paste0(
            "SNP重合率过低 (", round(snp_overlap_rate * 100, 2), "%)")
    } else {
        pgs_results$WARNING <- ""
    }
    
    # 添加风险等级评估
    pgs_results$RISK_LEVEL <- sapply(pgs_results$PERCENTILE, function(p) {
        if (is.na(p)) {
            return("无法评估")
        } else if (p >= 0.8) {
            return("高风险")
        } else if (p >= 0.6) {
            return("中高风险")
        } else if (p >= 0.4) {
            return("中等风险")
        } else if (p >= 0.2) {
            return("中低风险")
        } else {
            return("低风险")
        }
    })
    
    cat("PGS报告生成完成\n")
    
    return(pgs_results)
}
