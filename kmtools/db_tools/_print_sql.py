from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import get_lexer_by_name
from pygments.styles import get_style_by_name


def print_sql(sql_string, lexer="sql", style="default"):
    """Print SQL query with colorful formatting.

    Note
    ----
    This function is not tested because we don't want to
    bring in IPython as a dependency.
    """
    from IPython.display import display, HTML

    lexer = get_lexer_by_name(lexer)
    formatter = HtmlFormatter(style=get_style_by_name(style), noclasses=True)
    return display(HTML(highlight(sql_string, lexer, formatter)))
