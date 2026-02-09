#!/usr/bin/env python3
"""
Skill Collection Validator
全面检查所有 skill 的可运行性
"""

import os
import sys
import subprocess
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import yaml

class SkillValidator:
    def __init__(self, base_path):
        self.base_path = Path(base_path)
        self.skills_dir = self.base_path / "scientific-skills"
        self.results = []
        
    def validate_all(self):
        """验证所有 skill"""
        skills = [d for d in self.skills_dir.iterdir() if d.is_dir()]
        total = len(skills)
        
        print(f"🔍 开始验证 {total} 个 skills...\n")
        
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = {executor.submit(self.validate_skill, skill): skill for skill in skills}
            
            for i, future in enumerate(as_completed(futures)):
                skill = futures[future]
                try:
                    result = future.result()
                    self.results.append(result)
                    progress = (i + 1) / total * 100
                    status = "✅" if result['valid'] else "❌"
                    print(f"{status} [{i+1}/{total}] {skill.name} - {result['message']}")
                except Exception as e:
                    self.results.append({
                        'skill': skill.name,
                        'valid': False,
                        'message': f'验证异常: {str(e)}'
                    })
                    print(f"❌ [{i+1}/{total}] {skill.name} - 验证异常: {e}")
        
        self.print_summary()
        return self.results
    
    def validate_skill(self, skill_path):
        """验证单个 skill"""
        result = {
            'skill': skill_path.name,
            'valid': True,
            'checks': {},
            'message': '通过'
        }
        
        # 1. 检查 SKILL.md 存在
        skill_md = skill_path / "SKILL.md"
        if not skill_md.exists():
            result['valid'] = False
            result['checks']['skill_md'] = False
            result['message'] = '缺少 SKILL.md'
            return result
        result['checks']['skill_md'] = True
        
        # 2. 解析 SKILL.md YAML frontmatter
        try:
            with open(skill_md, 'r') as f:
                content = f.read()
            
            # 提取 YAML frontmatter
            if content.startswith('---'):
                parts = content.split('---', 2)
                if len(parts) >= 2:
                    yaml_content = parts[1]
                    metadata = yaml.safe_load(yaml_content)
                    result['checks']['yaml_valid'] = True
                    result['metadata'] = {
                        'name': metadata.get('name'),
                        'category': metadata.get('category'),
                        'skill_type': metadata.get('skill_type'),
                        'status': metadata.get('status')
                    }
        except Exception as e:
            result['checks']['yaml_valid'] = False
            result['checks']['yaml_error'] = str(e)
        
        # 3. 检查 scripts 目录和 main.py
        scripts_dir = skill_path / "scripts"
        main_py = scripts_dir / "main.py"
        
        if not scripts_dir.exists():
            result['valid'] = False
            result['checks']['scripts_dir'] = False
            result['message'] = '缺少 scripts 目录'
            return result
        result['checks']['scripts_dir'] = True
        
        if not main_py.exists():
            result['valid'] = False
            result['checks']['main_py'] = False
            result['message'] = '缺少 scripts/main.py'
            return result
        result['checks']['main_py'] = True
        
        # 4. Python 语法检查
        try:
            with open(main_py, 'r') as f:
                code = f.read()
            compile(code, str(main_py), 'exec')
            result['checks']['syntax'] = True
        except SyntaxError as e:
            result['valid'] = False
            result['checks']['syntax'] = False
            result['checks']['syntax_error'] = f"行 {e.lineno}: {e.msg}"
            result['message'] = f'语法错误: 行 {e.lineno}'
        
        # 5. 检查 requirements.txt（可选但有更好）
        req_file = skill_path / "requirements.txt"
        result['checks']['has_requirements'] = req_file.exists()
        
        # 6. 检查 references 目录（建议有）
        ref_dir = skill_path / "references"
        result['checks']['has_references'] = ref_dir.exists()
        
        # 7. 尝试导入检查（如果依赖简单）
        if result['valid'] and result['checks'].get('syntax', False):
            try:
                # 尝试检查是否有明显的导入错误
                result['checks']['import_check'] = self._check_imports(main_py)
            except Exception as e:
                result['checks']['import_check'] = f'导入检查失败: {str(e)}'
        
        return result
    
    def _check_imports(self, main_py):
        """检查 Python 文件中的导入语句"""
        with open(main_py, 'r') as f:
            content = f.read()
        
        imports = []
        for line in content.split('\n'):
            line = line.strip()
            if line.startswith('import ') or line.startswith('from '):
                imports.append(line)
        
        # 统计标准库 vs 第三方库
        std_libs = {'os', 'sys', 'json', 're', 'math', 'datetime', 'typing', 'pathlib', 
                   'collections', 'itertools', 'functools', 'hashlib', 'random', 'string'}
        
        external = []
        for imp in imports:
            module = imp.split()[1].split('.')[0]
            if module not in std_libs and not module.startswith('.'):
                try:
                    __import__(module)
                except ImportError:
                    external.append(module)
        
        return {
            'total_imports': len(imports),
            'missing_external': external
        }
    
    def print_summary(self):
        """打印汇总报告"""
        print("\n" + "="*60)
        print("📊 验证报告汇总")
        print("="*60)
        
        total = len(self.results)
        valid = sum(1 for r in self.results if r['valid'])
        invalid = total - valid
        
        print(f"\n总计: {total} skills")
        print(f"✅ 通过: {valid} ({valid/total*100:.1f}%)")
        print(f"❌ 失败: {invalid} ({invalid/total*100:.1f}%)")
        
        # 分类统计
        categories = {}
        for r in self.results:
            cat = r.get('metadata', {}).get('category', 'Unknown')
            if cat not in categories:
                categories[cat] = {'total': 0, 'valid': 0}
            categories[cat]['total'] += 1
            if r['valid']:
                categories[cat]['valid'] += 1
        
        print(f"\n📁 按分类统计:")
        for cat, stats in sorted(categories.items()):
            rate = stats['valid'] / stats['total'] * 100
            print(f"  {cat}: {stats['valid']}/{stats['total']} 通过 ({rate:.0f}%)")
        
        # 失败的详情
        if invalid > 0:
            print(f"\n❌ 失败的 skills:")
            for r in self.results:
                if not r['valid']:
                    print(f"  - {r['skill']}: {r['message']}")
        
        # 缺失依赖的警告
        missing_reqs = [r for r in self.results if not r['checks'].get('has_requirements', True)]
        if missing_reqs:
            print(f"\n⚠️  缺少 requirements.txt ({len(missing_reqs)} 个):")
            for r in missing_reqs[:10]:  # 只显示前10个
                print(f"  - {r['skill']}")
            if len(missing_reqs) > 10:
                print(f"  ... 还有 {len(missing_reqs) - 10} 个")
        
        # 生成详细 JSON 报告
        report_path = self.base_path / "validation_report.json"
        with open(report_path, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\n📄 详细报告已保存到: {report_path}")

if __name__ == "__main__":
    base_path = "/Users/z04030865/skills-collection"
    validator = SkillValidator(base_path)
    results = validator.validate_all()
    
    # 如果有失败，退出码非0
    failed = sum(1 for r in results if not r['valid'])
    sys.exit(1 if failed > 0 else 0)
