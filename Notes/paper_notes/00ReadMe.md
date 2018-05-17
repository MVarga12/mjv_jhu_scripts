## This directory contains notes on papers read

### Command to create PDF from Markdown

- Use pandoc to create a PDF from markdown via LaTeX 
  - Command:

  ```
  pandoc --from markdown file.md --pdf-engine=pdflatex -o file.pdf
  ```

- If you need math, use the extension `tex_math_single_backslash`, which treats \\(...\\) as latex inline math and
  \\[...\\] as latex block math

  ``` pandoc --from markdown+tex_math_single_backslash file.md --pdf-engine=pdflatex -o file.pdf
  ```

  - Note that the `+extension_name` after markdown means include this extension. To exclude an extension, use
  `-extension_name`
- You can also use some select LaTeX packages, such as `geometry`.
  - To do so, add `-V package_name:package_option`.
  - So using `geometry`, you'd write `geometry:margin=1in`, or whatever you want it to be
- Can also use YAML headers, but I haven't looked into this yet
