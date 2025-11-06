"""Example of well-structured Python module following best practices.

This module demonstrates:
- PEP 8 compliance
- Proper type hints
- Comprehensive docstrings
- Error handling
- Modular design
- DRY principle
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Protocol

# Constants
DEFAULT_TIMEOUT = 30
MAX_RETRIES = 3
VALID_STATUSES = frozenset(['pending', 'processing', 'completed', 'failed'])

# Configure logging
logger = logging.getLogger(__name__)


# Custom Exceptions
class DataProcessingError(Exception):
    """Base exception for data processing errors."""
    pass


class ValidationError(DataProcessingError):
    """Raised when data validation fails."""
    pass


class ResourceNotFoundError(DataProcessingError):
    """Raised when required resource is not found."""
    pass


# Protocol for dependency injection
class DataStore(Protocol):
    """Protocol defining interface for data storage."""

    def save(self, key: str, value: dict) -> None:
        """Save data to storage."""
        ...

    def load(self, key: str) -> dict:
        """Load data from storage."""
        ...


# Dataclass for configuration
@dataclass
class ProcessorConfig:
    """Configuration for data processor.

    Attributes:
        input_dir: Directory containing input files
        output_dir: Directory for processed files
        timeout: Operation timeout in seconds
        validate: Whether to validate data before processing
        allowed_extensions: Set of allowed file extensions
    """
    input_dir: Path
    output_dir: Path
    timeout: int = DEFAULT_TIMEOUT
    validate: bool = True
    allowed_extensions: set[str] = field(default_factory=lambda: {'.txt', '.csv', '.json'})

    def __post_init__(self):
        """Validate configuration after initialization."""
        if not self.input_dir.exists():
            raise ValueError(f"Input directory does not exist: {self.input_dir}")

        if self.timeout <= 0:
            raise ValueError(f"Timeout must be positive, got: {self.timeout}")

        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)


# Main processing class
class DataProcessor:
    """Process data files with validation and error handling.

    This class provides a robust interface for processing data files
    with built-in validation, error handling, and logging.

    Attributes:
        config: Processor configuration
        stats: Processing statistics
    """

    def __init__(self, config: ProcessorConfig, store: DataStore | None = None):
        """Initialize data processor.

        Args:
            config: Processor configuration
            store: Optional data store for persistence
        """
        self.config = config
        self.store = store
        self.stats = {
            'processed': 0,
            'failed': 0,
            'skipped': 0
        }
        logger.info(f"DataProcessor initialized with config: {config}")

    def process_file(self, filepath: Path) -> dict[str, any]:
        """Process a single file.

        Args:
            filepath: Path to file to process

        Returns:
            Dictionary containing processing results with keys:
                - 'status': Processing status
                - 'lines': Number of lines processed
                - 'timestamp': Processing timestamp

        Raises:
            ResourceNotFoundError: If file doesn't exist
            ValidationError: If file validation fails
            DataProcessingError: If processing fails
        """
        if not filepath.exists():
            raise ResourceNotFoundError(f"File not found: {filepath}")

        if self.config.validate and not self._validate_file(filepath):
            raise ValidationError(f"File validation failed: {filepath}")

        try:
            content = self._read_file(filepath)
            processed_data = self._process_content(content)
            output_path = self._write_output(filepath, processed_data)

            result = {
                'status': 'completed',
                'lines': len(content.splitlines()),
                'timestamp': datetime.now().isoformat(),
                'output': str(output_path)
            }

            self.stats['processed'] += 1
            logger.info(f"Successfully processed {filepath}")

            return result

        except Exception as e:
            self.stats['failed'] += 1
            logger.error(f"Failed to process {filepath}: {e}", exc_info=True)
            raise DataProcessingError(f"Processing failed for {filepath}") from e

    def process_directory(self) -> list[dict[str, any]]:
        """Process all files in input directory.

        Returns:
            List of processing results for each file

        Example:
            >>> config = ProcessorConfig(Path("input"), Path("output"))
            >>> processor = DataProcessor(config)
            >>> results = processor.process_directory()
            >>> print(f"Processed {len(results)} files")
        """
        results = []

        for filepath in self._get_files():
            try:
                result = self.process_file(filepath)
                results.append(result)
            except DataProcessingError as e:
                logger.warning(f"Skipping file {filepath}: {e}")
                self.stats['skipped'] += 1
                continue

        logger.info(f"Processing complete. Stats: {self.stats}")
        return results

    def get_statistics(self) -> dict[str, int]:
        """Get processing statistics.

        Returns:
            Dictionary with processing counts
        """
        return self.stats.copy()

    def reset_statistics(self) -> None:
        """Reset all statistics to zero."""
        for key in self.stats:
            self.stats[key] = 0

    def _validate_file(self, filepath: Path) -> bool:
        """Validate file before processing.

        Args:
            filepath: Path to file

        Returns:
            True if file is valid
        """
        # Check file extension
        if filepath.suffix not in self.config.allowed_extensions:
            logger.warning(
                f"Invalid extension {filepath.suffix}, "
                f"allowed: {self.config.allowed_extensions}"
            )
            return False

        # Check file is not empty
        if filepath.stat().st_size == 0:
            logger.warning(f"File is empty: {filepath}")
            return False

        return True

    def _read_file(self, filepath: Path) -> str:
        """Read file content.

        Args:
            filepath: Path to file

        Returns:
            File content as string

        Raises:
            DataProcessingError: If reading fails
        """
        try:
            return filepath.read_text(encoding='utf-8')
        except UnicodeDecodeError as e:
            raise DataProcessingError(f"Failed to decode {filepath}") from e
        except IOError as e:
            raise DataProcessingError(f"Failed to read {filepath}") from e

    def _process_content(self, content: str) -> str:
        """Process file content.

        Args:
            content: Raw file content

        Returns:
            Processed content
        """
        # Remove empty lines
        lines = [line.strip() for line in content.splitlines() if line.strip()]

        # Add processing metadata
        header = f"# Processed at {datetime.now().isoformat()}\n"
        footer = f"\n# Total lines: {len(lines)}"

        return header + '\n'.join(lines) + footer

    def _write_output(self, original_path: Path, content: str) -> Path:
        """Write processed content to output file.

        Args:
            original_path: Original file path
            content: Processed content

        Returns:
            Path to output file
        """
        output_path = self.config.output_dir / f"processed_{original_path.name}"

        try:
            output_path.write_text(content, encoding='utf-8')
            logger.debug(f"Wrote output to {output_path}")
            return output_path
        except IOError as e:
            raise DataProcessingError(f"Failed to write output: {output_path}") from e

    def _get_files(self) -> list[Path]:
        """Get list of files to process.

        Returns:
            Sorted list of file paths
        """
        files = []
        for ext in self.config.allowed_extensions:
            files.extend(self.config.input_dir.glob(f"*{ext}"))

        return sorted(files)


# Utility functions
def create_processor_from_config_file(config_path: Path) -> DataProcessor:
    """Create processor from configuration file.

    Args:
        config_path: Path to TOML configuration file

    Returns:
        Configured DataProcessor instance

    Example:
        >>> processor = create_processor_from_config_file(Path("config.toml"))
        >>> results = processor.process_directory()
    """
    import tomllib

    with config_path.open('rb') as f:
        config_data = tomllib.load(f)

    config = ProcessorConfig(
        input_dir=Path(config_data['input_dir']),
        output_dir=Path(config_data['output_dir']),
        timeout=config_data.get('timeout', DEFAULT_TIMEOUT),
        validate=config_data.get('validate', True)
    )

    return DataProcessor(config)


def main() -> None:
    """Main entry point for command-line usage."""
    import argparse

    parser = argparse.ArgumentParser(description="Process data files")
    parser.add_argument('input_dir', type=Path, help="Input directory")
    parser.add_argument('output_dir', type=Path, help="Output directory")
    parser.add_argument(
        '--timeout',
        type=int,
        default=DEFAULT_TIMEOUT,
        help=f"Timeout in seconds (default: {DEFAULT_TIMEOUT})"
    )
    parser.add_argument(
        '--no-validate',
        action='store_true',
        help="Skip validation"
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help="Verbose output"
    )

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Create processor
    config = ProcessorConfig(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        timeout=args.timeout,
        validate=not args.no_validate
    )

    processor = DataProcessor(config)

    # Process files
    try:
        results = processor.process_directory()
        stats = processor.get_statistics()

        print(f"\nProcessing complete!")
        print(f"Processed: {stats['processed']}")
        print(f"Failed: {stats['failed']}")
        print(f"Skipped: {stats['skipped']}")

    except Exception as e:
        logger.exception("Processing failed")
        print(f"Error: {e}")
        return


if __name__ == '__main__':
    main()
