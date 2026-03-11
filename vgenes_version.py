APP_NAME = "VGenes"
VERSION = (1, 1, 0)
SCHEMA_VERSION = 1

__version__ = ".".join(str(part) for part in VERSION)


def format_app_title(context=None):
    title = f"{APP_NAME} v{__version__}"
    if context:
        title = f"{title} - {context}"
    return title
