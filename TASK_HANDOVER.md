# Skills Collection 重构任务 - Session 交接文档

## 📊 当前进度

### ✅ 已完成
- **SKILLS_INDEX.md**: 215 个 skills 已按新 category 分类并排序
- **Category 同步**: 所有 SKILL.md 的 category 字段已更新
- **GitHub Commit**: 所有更改已 push
- **标杆 Skill**: `automated-soap-note-generator` 已按 K-Dense 标准重写 (700+ 行)

### 📁 项目位置
```
/Users/z04030865/skills-collection/
├── SKILLS_INDEX.md (已更新)
├── scientific-skills/
│   └── automated-soap-note-generator/SKILL.md (标杆 - 已重写)
└── reference-template/ (K-Dense 参考模板)
```

---

## 🎯 下一步任务

### 任务 1: 创建 3 个标杆 Skills（验证模板）
**目标**: 再写 3 个不同类型的 skill，确认模板稳定

**需要重写的 Skills**:
1. `western-blot-quantifier` (Wet Lab 类)
2. `clinical-data-cleaner` (Data 类)
3. `single-cell-rnaseq-pipeline` (Bioinfo 类)

**参考模板**: `/Users/z04030865/skills-collection/scientific-skills/automated-soap-note-generator/SKILL.md`

**标准结构**:
```markdown
---
name: [skill-name]
description: [一句话描述 + 使用边界]
allowed-tools: [Read, Write, Bash, Edit]
license: MIT
metadata:
    skill-author: AIPOCH
---

# [Title]

## Overview
[一句话总结 + Key Capabilities bullet points]

## When to Use
[✅ 什么时候用 / ❌ 什么时候不用 / 上下游 skill 配合]

## Integration with Other Skills
[Upstream / Downstream / Complete Workflow]

## Core Capabilities
[6 个模块，每个包含：功能说明、代码示例、参数表、Best Practices、Common Issues]

## Complete Workflow Example
[CLI 命令 + Python API + 预期输出]

## Quality Checklist
[Pre-analysis / During / Post-analysis / Before Publication]

## Common Pitfalls
[分类列出错误做法 → 正确做法]

## Troubleshooting
[实际问题 + 解决方案]

## References
[references/ 目录下的文档]

## Scripts
[scripts/ 目录下的文件]
```

---

### 任务 2: 写自动化脚本（批量升级）
**目标**: 自动化升级剩余 211 个 skills

**脚本功能**:
1. 读取现有 SKILL.md
2. 提取现有内容（description, usage, parameters）
3. 自动生成：
   - `## When to Use`（根据 category 和 name 推断）
   - `## Integration with Other Skills`（匹配上下游）
   - `## Quality Checklist`（标准模板）
   - `## Common Pitfalls`（根据 skill_type 匹配）
   - `## Troubleshooting`（常见问题模板）
4. 保留原有的代码示例和参数
5. 输出新的 SKILL.md

**技术方案**:
- 使用 Python + YAML 解析 frontmatter
- 使用 Jinja2 模板引擎
- 建立 skill 关系映射表（哪些 skill 是上下游）

---

### 任务 3: 批量生成并 Review
**目标**: 生成全部 215 个 skills 的新文档

**执行步骤**:
1. 跑自动化脚本生成初稿
2. 人工抽查 10-20 个关键 skills
3. 修复脚本问题
4. 批量 commit

---

## 📋 执行建议

### 方案 A: 先做标杆（推荐）
1. 人工重写 3 个标杆 skills（2-3 小时）
2. 基于标杆总结规律，写自动化脚本（4-6 小时）
3. 脚本批量生成 + review（1-2 天）

### 方案 B: 直接自动化
1. 直接写通用升级脚本（跳过标杆验证）
2. 批量生成（半天）
3. 人工 review 和修复（2-3 天）

**推荐方案 A**（质量更可控）

---

## 🔍 关键参考

### K-Dense 参考模板
```
/Users/z04030865/reference-template/scientific-skills/
├── rdkit/SKILL.md (Cheminformatics)
├── scanpy/SKILL.md (Single-cell analysis)
├── diffdock/SKILL.md (Molecular docking)
└── ... (其他 140 个 skills)
```

### 我们的标杆 Skill
```
/Users/z04030865/skills-collection/scientific-skills/automated-soap-note-generator/SKILL.md
```

### Skill 分类统计
```
AI/Tech: 2
Bioinfo: 25
Business: 1
Career: 9
Chemistry: 1
Clinical: 20
Data: 9
Education: 11
General: 13
Grant: 7
Info: 7
Integrity: 1
Management: 1
Operations: 6
Pharma: 29
Present: 6
Research: 20
Safety: 2
Utility: 7
Visual: 17
Wet Lab: 12
Writing: 9
```

---

## ⚡ 快速开始命令

```bash
# 进入项目目录
cd /Users/z04030865/skills-collection

# 查看当前 git 状态
git status

# 查看最近的 commit
git log --oneline -5

# 查看标杆 skill
cat scientific-skills/automated-soap-note-generator/SKILL.md

# 查看 K-Dense 参考
cat reference-template/scientific-skills/rdkit/SKILL.md
```

---

## ❓ 常见问题

**Q: 需要全部重写吗？**
A: 不需要。自动化脚本会保留原有代码示例，只补充缺失的 section。

**Q: 质量标准是什么？**
A: 对标 `/reference-template/scientific-skills/` 里的 K-Dense skills。

**Q: 时间预估？**
A: 
- 3 个标杆：2-3 小时
- 自动化脚本：4-6 小时
- 批量生成+review：1-2 天

---

## 🎯 成功标准

- [ ] 3 个标杆 skills 达到 K-Dense 水准
- [ ] 自动化脚本能正确生成 80%+ 的 skills
- [ ] 全部 215 个 skills 都有完整的 When to Use / Integration / Quality Checklist
- [ ] 所有 commit 已 push 到 GitHub

---

**创建时间**: 2026-02-09
**创建者**: opencode
**任务状态**: 进行中
**下一个 Session 应该**: 从"任务 1"开始，先写 3 个标杆 skills
