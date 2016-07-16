"""IPython Tools."""
from pygments import highlight
from pygments.lexers import SqlLexer
from pygments.formatters import HtmlFormatter
from IPython.display import display, HTML


def print_sql(sql_query, style='default'):
    """Print SQL code with syntax highlighting.

    To view all available color styles, run::

        from pygments.styles import get_all_styles
        print(sorted(get_all_styles()))
    """
    formatter = HtmlFormatter(style=style)
    display(HTML('<style type="text/css">{}</style>{}'.format(
        formatter.get_style_defs('.highlight'),
        highlight(sql_query, SqlLexer(), formatter))))
