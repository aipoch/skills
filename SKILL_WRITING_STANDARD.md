# SKILL.md 写作规范 v1.0
## AIPOCH Skills Collection 标准

---

## 一、Frontmatter 规范

### 必需字段
```yaml
---
name: [skill-name]
description: [一句话核心功能 + 使用边界说明]
allowed-tools: [Read, Write, Bash, Edit]  # 根据实际需要选择
license: MIT
metadata:
    skill-author: AIPOCH
---
```

### 禁止字段（已移除）
- ❌ version
- ❌ category（已在目录结构中体现）
- ❌ tags
- ❌ author
- ❌ status
- ❌ risk_level
- ❌ skill_type
- ❌ owner
- ❌ reviewer
- ❌ last_updated

---

## 二、文档结构标准（必须按顺序）

### 1. 标题（H1）
```markdown
# [Skill Name]
```

### 2. Overview（H2）
**内容要求：**
- 一句话总结核心功能
- 3-6 个 Key Capabilities（bullet points）

**格式：**
```markdown
## Overview

[一句话总结]

**Key Capabilities:**
- **功能1**: 简短说明
- **功能2**: 简短说明
- ...
```

### 3. When to Use（H2）
**必须包含三个部分：**

```markdown
## When to Use

**✅ Use this skill when:**
- 具体场景1
- 具体场景2
- ...

**❌ Do NOT use when:**
- 不适用场景1 → [替代skill]
- 不适用场景2 → [替代skill]
- ...

**Related Skills:**
- **上游**: [skill1], [skill2]
- **下游**: [skill3], [skill4]
```

### 4. Integration with Other Skills（H2）
**必须包含：**

```markdown
## Integration with Other Skills

**Upstream Skills:**
- `[skill-name]`: 上游功能说明

**Downstream Skills:**
- `[skill-name]`: 下游功能说明

**Complete Workflow:**
```
上游Skill → 本Skill → 下游Skill
```
```

### 5. Core Capabilities（H2）
**每个能力一个 H3 小节，包含：**

```markdown
### 1. [能力名称]

[功能说明]

```python
# 代码示例
from scripts.main import SkillName
skill = SkillName()
result = skill.function()
```

**Parameters:**
| Parameter | Type | Required | Description | Default |
|-----------|------|----------|-------------|---------|
| param1 | str | Yes | 说明 | None |

**Best Practices:**
- ✅ 应该这样做
- ✅ 另一个最佳实践

**Common Issues and Solutions:**

**Issue: 问题描述**
- Symptom: 症状
- Solution: 解决方案
```

**数量要求：**
- 最少 3 个，最多 8 个核心能力
- 每个能力必须包含：代码示例、参数表、Best Practices、Common Issues

### 6. Complete Workflow Example（H2）
**必须包含：**

```markdown
## Complete Workflow Example

**从输入到输出的完整流程：**

```bash
# Step 1: 步骤1
python scripts/main.py --param1 value1

# Step 2: 步骤2
python scripts/main.py --param2 value2
```

**Python API Usage:**

```python
from scripts.main import SkillName
from scripts.utils import Helper

# 完整代码示例
skill = SkillName()
result = skill.process()
```

**Expected Output Files:**
```
output/
├── file1.ext
└── file2.ext
```
```

### 7. Quality Checklist（H2）
**必须分阶段：**

```markdown
## Quality Checklist

**Pre-analysis Checks:**
- [ ] 检查项1
- [ ] 检查项2

**During [Process]:**
- [ ] 检查项1
- [ ] 检查项2

**Post-analysis Verification:**
- [ ] 检查项1
- [ ] 检查项2

**Before [Final Use]:**
- [ ] 检查项1
- [ ] 检查项2
```

**要求：**
- 至少 4 个阶段
- 每个阶段至少 3 个检查项
- 关键检查项标记 **CRITICAL**

### 8. Common Pitfalls（H2）
**分类列出，每类包含：**

```markdown
## Common Pitfalls

**[类别1] Issues:**
- ❌ **错误做法** → 后果
  - ✅ 正确做法

**[类别2] Issues:**
- ❌ **错误做法** → 后果
  - ✅ 正确做法
```

**要求：**
- 至少 3 个类别（如 Input Issues, Analysis Issues, Output Issues）
- 每个类别至少 3 个常见错误
- 必须提供 ✅ 正确做法

### 9. Troubleshooting（H2）
**问题-原因-解决方案格式：**

```markdown
## Troubleshooting

**Problem: [问题名称]**
- Symptoms: 症状描述
- Causes: 可能原因
- Solutions:
  - 解决方案1
  - 解决方案2
```

**要求：**
- 至少 5 个常见问题
- 每个问题必须包含 Symptoms, Causes, Solutions

### 10. References（H2）
```markdown
## References

Available in `references/` directory:

- `file1.md` - 说明
- `file2.md` - 说明
```

### 11. Scripts（H2）
```markdown
## Scripts

Located in `scripts/` directory:

- `main.py` - 主功能
- `utils.py` - 工具函数
```

### 12. 可选章节（根据需要添加）
- Performance and Resources
- Limitations
- Version History

---

## 三、写作质量标准

### 1. 语言规范
- **全英文**（包括代码注释）
- **专业术语准确**（医学/生物/化学术语）
- **简洁明了**，避免冗余

### 2. 格式规范
- 使用 `##` 和 `###` 正确分级
- 代码块使用 ```python 和 ```bash 标记
- 表格使用标准 Markdown 表格格式
- 重要提示使用 **粗体** 或 ⚠️ 符号

### 3. 代码示例规范
- 每个 Core Capability 必须有代码示例
- 代码必须可运行（语法正确）
- 包含必要的注释
- 展示常见用法和边界情况

### 4. 安全检查
- **⚠️ CRITICAL**: 必须强调安全警告（如医疗安全、数据安全）
- 所有医疗相关 skill 必须强调"需专业人员审核"
- 数据处理 skill 必须提及隐私保护

---

## 四、写作流程

### 每写 5 个 skills 必须休息
防止大脑疲劳导致质量下降

### 质量自查清单（写完每个 skill 后检查）

**Frontmatter:**
- [ ] name 正确
- [ ] description 包含使用边界
- [ ] allowed-tools 合理
- [ ] license 是 MIT
- [ ] metadata.skill-author 是 AIPOCH
- [ ] **没有多余的字段**

**内容结构:**
- [ ] Overview 有 Key Capabilities
- [ ] When to Use 有 ✅ ❌ 和 Related Skills
- [ ] Integration 有 Upstream/Downstream/Workflow
- [ ] Core Capabilities 至少 3 个，每个有代码+参数+Best Practices+Issues
- [ ] Complete Workflow 有 CLI + Python + Expected Output
- [ ] Quality Checklist 有 4+ 阶段，每阶段 3+ 项
- [ ] Common Pitfalls 有 3+ 类别，每类 3+ 错误
- [ ] Troubleshooting 有 5+ 问题
- [ ] References 和 Scripts 章节存在

**质量检查:**
- [ ] 全英文
- [ ] 代码语法正确
- [ ] 表格格式正确
- [ ] 没有中文标点符号（，。！）
- [ ] 链接和路径正确

---

## 五、特殊类别写作提示

### Clinical 类技能
- 必须强调医疗安全
- 必须说明需医生审核
- 必须提及 HIPAA 合规

### Wet Lab 类技能
- 必须包含实验安全提示
- 必须提供替代方案（如果实验失败）
- 必须说明仪器要求

### Data 类技能
- 必须说明输入数据格式要求
- 必须提及数据隐私保护
- 必须包含数据验证步骤

### Bioinfo 类技能
- 必须说明计算资源要求
- 必须提供参考文献
- 必须说明软件依赖版本

---

## 六、禁止事项

❌ **禁止：**
1. 使用中文（包括注释）
2. 添加 version/status/risk_level 等已移除字段
3. 跳过任何必需章节
4. 代码示例不可运行
5. Common Pitfalls 只列错误不给正确做法
6. 使用口语化表达
7. 复制粘贴不修改（必须定制化每个 skill）
8. 连续写超过 5 个不休息

---

## 七、示例参考

**标杆 Skill（必须严格对标）：**
`/Users/z04030865/skills-collection/scientific-skills/automated-soap-note-generator/SKILL.md`

**K-Dense 参考：**
`/Users/z04030865/reference-template/scientific-skills/rdkit/SKILL.md`

---

**规范版本**: v1.0
**创建日期**: 2026-02-09
**适用范围**: 所有 215 个 skills
**执行要求**: 必须严格遵守
