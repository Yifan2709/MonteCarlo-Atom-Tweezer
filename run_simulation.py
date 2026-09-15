"""仓库根目录入口；各层级仍使用各自的物理模型和配置。"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent /
                       "simulation/level2_joint_transfer/src"))
from simulation_workflow.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
