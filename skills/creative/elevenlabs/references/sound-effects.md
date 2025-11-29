# ElevenLabs Sound Effects Generation Reference

## Overview

ElevenLabs Sound Effects API generates high-quality sound effects from text descriptions. Perfect for game development, video production, UI/UX design, and content creation, the API creates realistic audio from simple prompts without requiring audio libraries or sample collections.

## Model

**Eleven Text-to-Sound v2**

**Model ID:** `eleven_text_to_sound_v2`

**Capabilities:**
- Text-to-sound-effect generation
- Customizable duration
- Prompt influence control
- High-quality audio output
- Suitable for professional use

## API Endpoint

**URL:** `POST https://api.elevenlabs.io/v1/sound-generation`

**Authentication:** API key in `xi-api-key` header or via SDK

## API Usage

### Basic Setup

```python
import os
from elevenlabs.client import ElevenLabs
from pathlib import Path

# Initialize client
client = ElevenLabs(api_key=os.environ.get("ELEVENLABS_API_KEY"))
```

### Simple Sound Effect

```python
# Generate a basic sound effect
audio = client.text_to_sound_effects.convert(
    text="footsteps on wooden floor, slow walking pace",
    duration_seconds=3.0
)

# Save to file
output_path = Path("footsteps.mp3")
with output_path.open("wb") as f:
    for chunk in audio:
        f.write(chunk)

print(f"Sound effect saved: {output_path}")
```

### With All Parameters

```python
# Generate with full parameter control
audio = client.text_to_sound_effects.convert(
    text="large explosion with debris falling, action movie style",
    duration_seconds=5.0,
    prompt_influence=0.5  # How closely to follow the prompt
)

output_path = Path("explosion.mp3")
with output_path.open("wb") as f:
    for chunk in audio:
        f.write(chunk)
```

## Parameters

### Required Parameters

**text** (string)
- Description of the desired sound effect
- Be specific and descriptive
- Include context, mood, and characteristics
- Examples:
  - "gentle rain falling on leaves"
  - "car engine starting, diesel truck"
  - "sword clashing against metal shield"

### Optional Parameters

**duration_seconds** (number)
- Length of the generated audio in seconds
- Typical range: 0.5 to 22 seconds
- Default varies by model
- Longer durations take more time to generate

**prompt_influence** (number, 0.0 - 1.0)
- Controls how literally the model interprets the prompt
- **0.0 - 0.3**: More creative interpretation, varied results
- **0.4 - 0.6**: Balanced creativity and adherence (recommended)
- **0.7 - 1.0**: Strict adherence to prompt, literal interpretation
- Default: 0.5

**output_format** (string)
- Audio codec and quality
- Format: `{codec}_{sample_rate}_{bitrate}`
- Examples: `mp3_44100_128`, `mp3_22050_64`
- Default: `mp3_44100_128`

## Use Cases

### Game Audio

**Combat Sounds:**
```python
combat_sounds = [
    ("sword_swing", "sword whooshing through air, fast swing"),
    ("shield_block", "metal shield blocking attack, metallic clang"),
    ("arrow_hit", "arrow hitting wooden target with thud"),
    ("magic_spell", "magical energy spell casting, ethereal whoosh")
]

for name, description in combat_sounds:
    audio = client.text_to_sound_effects.convert(
        text=description,
        duration_seconds=1.5
    )
    save_audio(audio, f"game_{name}.mp3")
```

**Environmental Sounds:**
```python
environment = [
    ("forest_ambience", "forest birds chirping, leaves rustling, nature sounds"),
    ("dungeon_drip", "water dripping in dark stone dungeon, echo"),
    ("fire_crackle", "campfire crackling and burning, warm flames"),
    ("wind_howl", "strong wind howling through cave, eerie atmosphere")
]

for name, description in environment:
    audio = client.text_to_sound_effects.convert(
        text=description,
        duration_seconds=5.0,
        prompt_influence=0.4  # Allow creative interpretation
    )
    save_audio(audio, f"env_{name}.mp3")
```

**UI Sounds:**
```python
ui_sounds = [
    ("button_click", "soft UI button click, pleasant tone", 0.3),
    ("success_chime", "success notification sound, positive chime", 0.8),
    ("error_buzz", "error notification, warning buzz", 0.5),
    ("menu_open", "menu sliding open, smooth swoosh", 0.6)
]

for name, description, duration in ui_sounds:
    audio = client.text_to_sound_effects.convert(
        text=description,
        duration_seconds=duration
    )
    save_audio(audio, f"ui_{name}.mp3")
```

### Video Production

**Transitions:**
```python
transitions = {
    "whoosh": "fast whoosh transition sound, energetic swipe",
    "glitch": "digital glitch transition, electronic distortion",
    "rise": "rising tension sound, building anticipation",
    "impact": "dramatic impact hit, powerful punch"
}

for name, description in transitions.items():
    audio = client.text_to_sound_effects.convert(
        text=description,
        duration_seconds=1.0
    )
    save_audio(audio, f"transition_{name}.mp3")
```

**Ambiences:**
```python
# Longer ambient sounds for backgrounds
ambiences = [
    ("city_street", "busy city street ambience, traffic and pedestrians", 10.0),
    ("office", "office ambience, keyboard typing and phone ringing", 10.0),
    ("cafe", "coffee shop ambience, chatter and espresso machine", 10.0),
    ("nature", "peaceful nature ambience, birds and gentle breeze", 10.0)
]

for name, description, duration in ambiences:
    audio = client.text_to_sound_effects.convert(
        text=description,
        duration_seconds=duration
    )
    save_audio(audio, f"ambience_{name}.mp3")
```

### Podcast/Audio Production

**Production Elements:**
```python
production = {
    "intro_swoosh": "podcast intro whoosh, professional and clean",
    "transition_pop": "section transition pop, upbeat and modern",
    "outro_fade": "podcast outro background fade, smooth ending",
    "notification": "listener notification sound, attention-grabbing ding"
}

for name, description in production.items():
    audio = client.text_to_sound_effects.convert(
        text=description,
        duration_seconds=2.0
    )
    save_audio(audio, f"podcast_{name}.mp3")
```

## Prompt Engineering for Sound Effects

### Effective Prompt Structure

**Formula:** `[Sound Source] + [Action/Movement] + [Context/Environment] + [Mood/Style]`

**Examples:**

**Basic:**
```
"door creaking open"
```

**Better:**
```
"old wooden door creaking open slowly"
```

**Best:**
```
"old wooden door creaking open slowly in abandoned house, horror atmosphere, eerie"
```

### Prompt Components

**1. Sound Source**
- Be specific about the object/source
- Examples: "wooden door", "metal gate", "glass window"

**2. Action/Movement**
- Describe what's happening
- Examples: "creaking open", "slamming shut", "shattering"

**3. Context/Environment**
- Where the sound occurs
- Examples: "in small room", "outdoor echo", "underwater"

**4. Characteristics**
- Speed: slow, fast, rapid
- Intensity: gentle, loud, powerful
- Quality: crisp, muffled, distorted

**5. Mood/Atmosphere**
- Emotional context
- Examples: "scary", "cheerful", "dramatic", "peaceful"

**6. Style/Genre**
- Production style
- Examples: "action movie", "horror film", "cartoon", "realistic"

### Examples by Category

**Footsteps:**
```python
footsteps = [
    "footsteps on gravel path, slow walking pace, outdoor",
    "high heels clicking on marble floor, confident stride, indoor echo",
    "heavy boots on metal grating, industrial setting",
    "bare feet on wet sand, beach ambience"
]
```

**Doors:**
```python
doors = [
    "car door closing firmly, modern vehicle",
    "elevator doors sliding open with electronic chime",
    "jail cell door slamming shut, metallic clang, dramatic",
    "automatic sliding doors opening, futuristic sci-fi"
]
```

**Nature:**
```python
nature = [
    "thunder rumbling in distance, storm approaching",
    "ocean waves crashing on rocky shore, powerful",
    "gentle stream flowing over rocks, peaceful forest",
    "strong wind through pine trees, mountain atmosphere"
]
```

**Mechanical:**
```python
mechanical = [
    "car engine revving, sports car, powerful acceleration",
    "factory machinery operating, industrial metallic sounds",
    "clock ticking, antique grandfather clock, steady rhythm",
    "helicopter rotors spinning, close flyby, dramatic"
]
```

**Fantasy/Sci-Fi:**
```python
fantasy_scifi = [
    "magical portal opening, ethereal energy swirl, mystical",
    "laser gun firing, sci-fi weapon, futuristic",
    "dragon roaring, powerful and intimidating, fantasy",
    "spaceship engine humming, deep bass, sci-fi ambience"
]
```

### Advanced Prompting Techniques

**Layered Descriptions:**
```python
# Multiple sound elements in one prompt
audio = client.text_to_sound_effects.convert(
    text="thunderstorm with rain falling, wind blowing, and distant thunder rumbles",
    duration_seconds=8.0
)
```

**Temporal Progression:**
```python
# Describe sound evolution over time
audio = client.text_to_sound_effects.convert(
    text="footsteps approaching from distance, getting louder and closer, then stopping",
    duration_seconds=5.0
)
```

**Comparative References:**
```python
# Reference familiar sounds
audio = client.text_to_sound_effects.convert(
    text="explosion similar to action movie, debris falling, dramatic impact",
    duration_seconds=4.0
)
```

## Duration Guidelines

### Short Sounds (0.3 - 2 seconds)
- UI clicks and feedback
- Button presses
- Notifications
- Single impacts
- Quick transitions

```python
short_sfx = [
    ("click", "UI button click", 0.3),
    ("ding", "notification bell", 0.5),
    ("pop", "bubble pop sound", 0.4),
    ("snap", "finger snap", 0.3)
]
```

### Medium Sounds (2 - 5 seconds)
- Footsteps sequences
- Door actions
- Vehicle sounds
- Object interactions
- Scene transitions

```python
medium_sfx = [
    ("footsteps", "footsteps on wooden floor, 3 steps", 2.0),
    ("door", "door opening and closing", 3.0),
    ("car_start", "car engine starting up", 4.0)
]
```

### Long Sounds (5 - 22 seconds)
- Ambient loops
- Environmental sounds
- Continuous processes
- Atmospheric backgrounds

```python
long_sfx = [
    ("rain", "gentle rain falling, ambient loop", 10.0),
    ("crowd", "crowded marketplace ambience", 15.0),
    ("machinery", "factory machinery operating continuously", 20.0)
]
```

## Output Formats

### Recommended Formats

**High Quality (Default):**
```python
audio = client.text_to_sound_effects.convert(
    text="explosion",
    duration_seconds=3.0,
    output_format="mp3_44100_128"  # 44.1kHz, 128kbps
)
```

**Standard Quality:**
```python
audio = client.text_to_sound_effects.convert(
    text="footsteps",
    duration_seconds=2.0,
    output_format="mp3_44100_96"  # 44.1kHz, 96kbps
)
```

**Compact (Smaller Files):**
```python
audio = client.text_to_sound_effects.convert(
    text="ui click",
    duration_seconds=0.5,
    output_format="mp3_22050_64"  # 22.05kHz, 64kbps
)
```

## Batch Generation

### Generate Multiple Effects

```python
def batch_generate_sfx(effects_list):
    """Generate multiple sound effects from a list"""
    results = []

    for name, description, duration in effects_list:
        print(f"Generating: {name}...")

        try:
            audio = client.text_to_sound_effects.convert(
                text=description,
                duration_seconds=duration
            )

            output_path = Path(f"sfx_{name}.mp3")
            with output_path.open("wb") as f:
                for chunk in audio:
                    f.write(chunk)

            results.append((name, output_path, "success"))
            print(f"  ✓ Saved: {output_path}")

        except Exception as e:
            results.append((name, None, f"error: {e}"))
            print(f"  ✗ Failed: {e}")

    return results

# Use the batch generator
game_sfx = [
    ("jump", "character jumping, game sound effect", 1.0),
    ("land", "character landing on ground, impact", 0.8),
    ("collect", "collecting coin or item, positive chime", 0.6),
    ("damage", "taking damage, hurt sound", 0.7)
]

results = batch_generate_sfx(game_sfx)
```

### Organized Export

```python
from pathlib import Path

def organize_sfx_library(sfx_categories):
    """Generate and organize sound effects by category"""

    base_dir = Path("sfx_library")
    base_dir.mkdir(exist_ok=True)

    for category, effects in sfx_categories.items():
        category_dir = base_dir / category
        category_dir.mkdir(exist_ok=True)

        for name, description, duration in effects:
            audio = client.text_to_sound_effects.convert(
                text=description,
                duration_seconds=duration
            )

            output_path = category_dir / f"{name}.mp3"
            with output_path.open("wb") as f:
                for chunk in audio:
                    f.write(chunk)

            print(f"Created: {category}/{name}.mp3")

# Organize by category
sfx_library = {
    "ui": [
        ("button_click", "UI button click", 0.3),
        ("hover", "UI hover sound", 0.2),
    ],
    "environment": [
        ("wind", "gentle wind blowing", 5.0),
        ("rain", "rain falling", 8.0),
    ],
    "combat": [
        ("sword_hit", "sword hitting armor", 0.8),
        ("shield_block", "blocking with shield", 0.7),
    ]
}

organize_sfx_library(sfx_library)
```

## Best Practices

### 1. Descriptive Prompts
- Be specific and detailed
- Include context and atmosphere
- Mention technical characteristics
- Reference style or genre when helpful

### 2. Appropriate Duration
- Match duration to sound type
- UI sounds: 0.3-0.8 seconds
- Actions: 1-3 seconds
- Ambiences: 5-20 seconds
- Consider loop length for repeating sounds

### 3. Prompt Influence
- Start with 0.5 (balanced)
- Increase for more literal interpretation
- Decrease for more creative variation
- Experiment to find sweet spot

### 4. Iteration
- Generate multiple variations
- Try different descriptions
- Adjust prompt_influence
- Compare results

### 5. Organization
- Use clear naming conventions
- Organize by category or use case
- Document prompts used
- Track which settings worked best

### 6. Quality Control
- Listen to generated audio
- Check duration matches needs
- Verify audio quality is sufficient
- Test in target application

## Common Issues & Solutions

### Issue: Sound doesn't match description

**Solutions:**
- Make prompt more specific
- Increase prompt_influence
- Add more descriptive details
- Try different phrasing

### Issue: Too short/long

**Solutions:**
- Adjust duration_seconds parameter
- For very short sounds, use 0.3-0.5 seconds minimum
- For loops, generate longer and trim if needed

### Issue: Poor quality

**Solutions:**
- Use higher bitrate format
- Increase sample rate
- Generate longer duration
- Check source prompt clarity

### Issue: Inconsistent results

**Solutions:**
- Use higher prompt_influence for consistency
- Make prompts more specific
- Generate multiple samples and choose best
- Use same settings for related sounds

## Advanced Techniques

### Creating Sound Variations

```python
# Generate variations of the same sound
base_prompt = "footsteps on wooden floor"

variations = [
    f"{base_prompt}, slow walking",
    f"{base_prompt}, fast running",
    f"{base_prompt}, sneaking quietly",
    f"{base_prompt}, heavy stomping"
]

for i, prompt in enumerate(variations):
    audio = client.text_to_sound_effects.convert(
        text=prompt,
        duration_seconds=2.0
    )
    save_audio(audio, f"footsteps_var{i+1}.mp3")
```

### Layering for Complexity

```python
# Generate separate layers to mix later
layers = [
    ("base", "deep bass rumble, low frequency", 5.0),
    ("mid", "metallic clanging and impacts", 5.0),
    ("high", "debris falling and scattering", 5.0)
]

for layer_name, description, duration in layers:
    audio = client.text_to_sound_effects.convert(
        text=description,
        duration_seconds=duration
    )
    save_audio(audio, f"explosion_{layer_name}.mp3")

# Mix layers using audio editing software or pydub
```

### Testing Different Influences

```python
# Test different prompt_influence values
prompt = "sword swinging through air"
influences = [0.3, 0.5, 0.7]

for influence in influences:
    audio = client.text_to_sound_effects.convert(
        text=prompt,
        duration_seconds=1.5,
        prompt_influence=influence
    )
    save_audio(audio, f"sword_swing_{int(influence*10)}.mp3")
```

## API Limits

### Rate Limits
- Subject to account tier limits
- Free tier: Limited requests per month
- Paid tier: Higher limits based on plan

### Character Limits
- Sound effect descriptions count toward character quota
- Longer prompts use more characters
- Monitor usage via dashboard

### Duration Limits
- Maximum duration varies by tier
- Typically up to 22 seconds
- Longer sounds may require paid tier

## Resources

- **API Reference**: https://elevenlabs.io/docs/api-reference/text-to-sound-effects
- **Sound Effects Guide**: https://elevenlabs.io/docs/capabilities/sound-effects
- **Usage Dashboard**: https://elevenlabs.io/app/usage
- **Python SDK**: https://github.com/elevenlabs/elevenlabs-python
- **Community Examples**: https://elevenlabs.io/showcase
