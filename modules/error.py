

class Error(Exception):
    """Base class for other exceptions"""
    pass


class DependencyNotFoundError(Error):
    """Raised when a dependency is not found."""
    def __init__(self, dependency_name, message="Dependency not found,"
                                                "check gdock.ini:"):
        self.dependency_name = dependency_name
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message} {self.dependency_name}'


class DependencyNotDefinedError(Error):
    """Raised when a needed third-party dependency is not properly defined."""
    def __init__(self, dependency_name, message="Dependency not defined in"
                                                "gdock.ini:"):
        self.dependency_name = dependency_name
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message} {self.dependency_name}'


class SectionNotDefinedError(Error):
    """Raised when a needed section is not properly defined."""
    def __init__(self, section_name, message="Section not defined in"
                                             "gdock.ini:"):
        self.section_name = section_name
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message} {self.section_name}'
